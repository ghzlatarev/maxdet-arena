"use client";

import { useCallback, useEffect, useState } from "react";

type EthereumProvider = {
  request: (request: { method: string; params?: unknown[] }) => Promise<unknown>;
};

declare global {
  interface Window {
    ethereum?: EthereumProvider;
  }
}

export type MaxDetBountyDeployment = {
  network: string;
  chain_id: number;
  status: string;
  contract_address: string | null;
  deployment_transaction: string | null;
  asset: string;
  minimum_winning_determinant: string;
  expected_runtime_codehash: string;
  source: string;
  explorer_url: string;
  rpc_url: string;
};

type ChainSnapshot = {
  balance: bigint;
  solved: boolean;
  deferredPrize: bigint;
};

const solvedSelector = "0x799320bb";
const minimumSelector = "0xb531b874";
const orderSelector = "0x4d5e07fb";
const claimablePrizeSelector = "0x22629a99";

function isAddress(value: string | null): value is string {
  return typeof value === "string" && /^0x[0-9a-fA-F]{40}$/.test(value);
}

function isExpectedChain(value: unknown, chainId: number): boolean {
  if (typeof value !== "string") return false;
  try {
    return BigInt(value) === BigInt(chainId);
  } catch {
    return false;
  }
}

function providerErrorCode(error: unknown): number | undefined {
  if (typeof error !== "object" || error === null || !("code" in error)) {
    return undefined;
  }
  const code = (error as { code?: unknown }).code;
  return typeof code === "number" ? code : undefined;
}

function parseEther(value: string): bigint {
  if (!/^\d+(?:\.\d{0,18})?$/.test(value)) {
    throw new Error("Enter a valid ETH amount with at most 18 decimals.");
  }
  const [whole, fraction = ""] = value.split(".");
  const wei = BigInt(whole) * 10n ** 18n;
  return wei + BigInt((fraction + "0".repeat(18)).slice(0, 18));
}

function formatEther(value: bigint): string {
  const whole = value / 10n ** 18n;
  const fraction = (value % 10n ** 18n)
    .toString()
    .padStart(18, "0")
    .slice(0, 4)
    .replace(/0+$/, "");
  return fraction ? `${whole}.${fraction}` : whole.toString();
}

function shortAddress(address: string): string {
  return `${address.slice(0, 6)}…${address.slice(-4)}`;
}

async function rpc<T>(
  rpcUrl: string,
  method: string,
  params: unknown[],
): Promise<T> {
  const response = await fetch(rpcUrl, {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: JSON.stringify({ jsonrpc: "2.0", id: 1, method, params }),
  });
  if (!response.ok) throw new Error(`Sepolia RPC returned ${response.status}`);
  const payload = (await response.json()) as {
    result?: T;
    error?: { message?: string };
  };
  if (payload.error || payload.result === undefined) {
    throw new Error(payload.error?.message ?? "Sepolia RPC request failed");
  }
  return payload.result;
}

export function MaxDetBountyCard({
  deployment,
  rulesUrl,
  sourceUrl,
}: {
  deployment: MaxDetBountyDeployment;
  rulesUrl: string;
  sourceUrl: string;
}) {
  const address =
    deployment.status === "live" && isAddress(deployment.contract_address)
      ? deployment.contract_address
      : null;
  const [snapshot, setSnapshot] = useState<ChainSnapshot | null>(null);
  const [amount, setAmount] = useState("0.01");
  const [message, setMessage] = useState("");
  const [verificationError, setVerificationError] = useState("");
  const [transactionHash, setTransactionHash] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const refresh = useCallback(async () => {
    if (!address) {
      setSnapshot(null);
      return;
    }
    setSnapshot(null);
    try {
      const [
        chainHex,
        codeHex,
        balanceHex,
        solvedHex,
        minimumHex,
        orderHex,
        deferredPrizeHex,
      ] = await Promise.all([
        rpc<string>(deployment.rpc_url, "eth_chainId", []),
        rpc<string>(deployment.rpc_url, "eth_getCode", [address, "latest"]),
        rpc<string>(deployment.rpc_url, "eth_getBalance", [address, "latest"]),
        rpc<string>(deployment.rpc_url, "eth_call", [
          { to: address, data: solvedSelector },
          "latest",
        ]),
        rpc<string>(deployment.rpc_url, "eth_call", [
          { to: address, data: minimumSelector },
          "latest",
        ]),
        rpc<string>(deployment.rpc_url, "eth_call", [
          { to: address, data: orderSelector },
          "latest",
        ]),
        rpc<string>(deployment.rpc_url, "eth_call", [
          { to: address, data: claimablePrizeSelector },
          "latest",
        ]),
      ]);

      if (!isExpectedChain(chainHex, deployment.chain_id)) {
        throw new Error("RPC is connected to the wrong chain");
      }
      if (/^0x0*$/.test(codeHex)) {
        throw new Error("no contract code exists at the published address");
      }
      const codehash = await rpc<string>(deployment.rpc_url, "web3_sha3", [
        codeHex,
      ]);
      if (
        codehash.toLowerCase() !==
        deployment.expected_runtime_codehash.toLowerCase()
      ) {
        throw new Error("on-chain bytecode does not match the release registry");
      }
      if (BigInt(minimumHex) !== BigInt(deployment.minimum_winning_determinant)) {
        throw new Error("on-chain winning threshold does not match the registry");
      }
      if (BigInt(orderHex) !== 23n) {
        throw new Error("on-chain matrix order is not 23");
      }

      setSnapshot({
        balance: BigInt(balanceHex),
        solved: BigInt(solvedHex) !== 0n,
        deferredPrize: BigInt(deferredPrizeHex),
      });
      setVerificationError("");
    } catch (error) {
      setSnapshot(null);
      setVerificationError(
        error instanceof Error
          ? `Contract check failed: ${error.message}.`
          : "Contract check failed.",
      );
    }
  }, [
    address,
    deployment.chain_id,
    deployment.expected_runtime_codehash,
    deployment.minimum_winning_determinant,
    deployment.rpc_url,
  ]);

  useEffect(() => {
    void refresh();
    const interval = window.setInterval(() => void refresh(), 30_000);
    return () => window.clearInterval(interval);
  }, [refresh]);

  async function donate() {
    if (!address || !snapshot || snapshot.solved) return;
    const provider = window.ethereum;
    if (!provider) {
      setMessage("Open this page in an Ethereum wallet to donate test ETH.");
      return;
    }
    setBusy(true);
    setMessage("");
    try {
      const value = parseEther(amount);
      if (value === 0n) throw new Error("Donation must be greater than zero.");

      const chainId = `0x${deployment.chain_id.toString(16)}`;
      const currentChain = await provider.request({
        method: "eth_chainId",
      });
      if (!isExpectedChain(currentChain, deployment.chain_id)) {
        try {
          await provider.request({
            method: "wallet_switchEthereumChain",
            params: [{ chainId }],
          });
        } catch (error) {
          if (providerErrorCode(error) !== 4902) throw error;
          await provider.request({
            method: "wallet_addEthereumChain",
            params: [
              {
                chainId,
                chainName: deployment.network,
                nativeCurrency: { name: "Sepolia ETH", symbol: "ETH", decimals: 18 },
                rpcUrls: [deployment.rpc_url],
                blockExplorerUrls: [deployment.explorer_url],
              },
            ],
          });
          await provider.request({
            method: "wallet_switchEthereumChain",
            params: [{ chainId }],
          });
        }
      }

      const accounts = (await provider.request({
        method: "eth_requestAccounts",
      })) as string[];
      const activeChain = await provider.request({ method: "eth_chainId" });
      if (!isExpectedChain(activeChain, deployment.chain_id)) {
        throw new Error("Wallet is not connected to Sepolia.");
      }
      if (!Array.isArray(accounts) || !isAddress(accounts[0] ?? null)) {
        throw new Error("No valid wallet account is available.");
      }
      const hash = (await provider.request({
        method: "eth_sendTransaction",
        params: [
          {
            from: accounts[0],
            to: address,
            value: `0x${value.toString(16)}`,
          },
        ],
      })) as string;
      setTransactionHash(hash);
      setMessage("Donation sent. The balance updates after it is mined.");
      window.setTimeout(() => void refresh(), 12_000);
    } catch (error) {
      setMessage(error instanceof Error ? error.message : "Donation cancelled.");
    } finally {
      setBusy(false);
    }
  }

  const live = Boolean(address && snapshot);
  const solved = snapshot?.solved ?? false;
  const payoutPending = solved && (snapshot?.deferredPrize ?? 0n) > 0n;
  const status =
    deployment.status !== "live"
      ? "Ready to deploy"
      : !address
        ? "Invalid deployment"
        : !snapshot
          ? "Checking"
          : payoutPending
            ? "Payout pending"
            : solved
              ? "Paid"
              : "Open";

  return (
    <section className="maxdet-bounty" id="bounty" aria-labelledby="bounty-title">
      <header>
        <div>
          <span className="maxdet-bounty-network">Ethereum · Sepolia</span>
          <h3 id="bounty-title">On-chain bounty</h3>
        </div>
        <span className={`maxdet-bounty-status ${solved ? "is-solved" : ""}`}>
          <i aria-hidden="true" /> {status}
        </span>
      </header>

      <div className="maxdet-bounty-balance">
        <span>{payoutPending ? "Deferred prize" : "Prize balance"}</span>
        <strong>
          {address
            ? snapshot
              ? `${formatEther(snapshot.balance)} ETH`
              : "Loading…"
            : "Awaiting Sepolia deployment"}
        </strong>
        <p>
          {payoutPending
            ? "The winning matrix is verified. Its winner can redirect the deferred payout."
            : "Anyone can donate. The first successful committed reveal above the fixed frontier earns the entire balance."}
        </p>
      </div>

      <dl className="maxdet-bounty-facts">
        <div>
          <dt>Winning score</dt>
          <dd>
            {BigInt(deployment.minimum_winning_determinant).toLocaleString(
              "en-US",
            )}
          </dd>
        </div>
        <div>
          <dt>Payout</dt>
          <dd>100% · one shot</dd>
        </div>
        <div>
          <dt>Verifier</dt>
          <dd>Exact · fully on-chain</dd>
        </div>
      </dl>

      <div className="maxdet-bounty-actions">
        <div className="maxdet-donation-input">
          <input
            aria-label="Sepolia ETH donation amount"
            disabled={!live || solved || busy}
            inputMode="decimal"
            onChange={(event) => setAmount(event.target.value)}
            value={amount}
          />
          <span>ETH</span>
        </div>
        <button disabled={!live || solved || busy} onClick={() => void donate()} type="button">
          {busy
            ? "Confirming…"
            : payoutPending
              ? "Payout pending"
              : solved
                ? "Bounty paid"
                : live
                  ? "Donate test ETH"
                  : deployment.status === "live"
                    ? "Checking contract…"
                    : "Awaiting deployment"}
        </button>
      </div>

      <footer className="maxdet-bounty-footer">
        <span>Testnet only · irrevocable · no refunds</span>
        <nav aria-label="Bounty links">
          <a href={rulesUrl}>Rules ↗</a>
          <a href={address ? `${deployment.explorer_url}/address/${address}` : sourceUrl}>
            {address ? shortAddress(address) : "Contract source"} ↗
          </a>
        </nav>
      </footer>

      {message ? <p className="maxdet-bounty-message">{message}</p> : null}
      {verificationError ? (
        <p className="maxdet-bounty-message">{verificationError}</p>
      ) : null}
      {transactionHash ? (
        <a
          className="maxdet-bounty-transaction"
          href={`${deployment.explorer_url}/tx/${transactionHash}`}
        >
          View donation transaction ↗
        </a>
      ) : null}
    </section>
  );
}
