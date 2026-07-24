import fs from "node:fs";
import path from "node:path";
import { CopyCommand } from "@/components/CopyCommand";
import { MatrixField } from "@/components/MatrixField";
import challenge from "@/challenge.json";
import frontier from "@/data/frontier.json";

const githubUrl = "https://github.com/ghzlatarev/maxdet-arena";
const command =
  'git clone https://github.com/ghzlatarev/maxdet-arena.git && cd maxdet-arena && codex "Read AGENTS.md. Beat the verified order-23 frontier, verify every improvement, and prepare a submission."';
const ehlichBoundSquared = BigInt(
  challenge.verification.order_specific_bound.determinant_squared,
);

function loadReferenceMatrix(): number[][] {
  const matrixPath = path.join(
    process.cwd(),
    "references",
    "orrick-et-al-2003",
    "matrix.txt",
  );
  return fs
    .readFileSync(matrixPath, "utf8")
    .trim()
    .split(/\r?\n/)
    .map((line) => line.trim().split(/\s+/).map(Number));
}

function formatInteger(value: string): string {
  return BigInt(value).toLocaleString("en-US");
}

function loadTargetToBeat(): {
  absoluteDeterminant: string;
  source: string;
  detail: string;
} {
  let best = {
    absoluteDeterminant: frontier.target_to_beat.absolute_determinant,
    source: "Orrick et al. · 2003",
    detail: "2²² × 3 × 5⁶ × 67 × 211",
  };
  if (
    BigInt(frontier.arena_best.absolute_determinant) >
    BigInt(best.absoluteDeterminant)
  ) {
    best = {
      absoluteDeterminant: frontier.arena_best.absolute_determinant,
      source: frontier.arena_best.label,
      detail: `receipt ${frontier.arena_best.receipt_sha256.slice(0, 12)}…`,
    };
  }

  const submissions = path.join(process.cwd(), "submissions");
  if (!fs.existsSync(submissions)) return best;
  for (const handle of fs.readdirSync(submissions, { withFileTypes: true })) {
    if (!handle.isDirectory()) continue;
    const handlePath = path.join(submissions, handle.name);
    for (const result of fs.readdirSync(handlePath, { withFileTypes: true })) {
      if (!result.isDirectory()) continue;
      const receipt = JSON.parse(
        fs.readFileSync(
          path.join(handlePath, result.name, "receipt.json"),
          "utf8",
        ),
      ) as {
        receipt_sha256: string;
        score: { absolute_determinant: string };
      };
      if (
        BigInt(receipt.score.absolute_determinant) >
        BigInt(best.absoluteDeterminant)
      ) {
        best = {
          absoluteDeterminant: receipt.score.absolute_determinant,
          source: `${handle.name}/${result.name}`,
          detail: `receipt ${receipt.receipt_sha256.slice(0, 12)}…`,
        };
      }
    }
  }
  return best;
}

export default function Home() {
  const matrix = loadReferenceMatrix();
  const targetToBeat = loadTargetToBeat();
  const target = targetToBeat.absoluteDeterminant;
  const arenaBest = frontier.arena_best.absolute_determinant;
  const statusLabel =
    frontier.status === "private-dogfooding"
      ? "Private dogfooding"
      : "Open arena";
  const progress =
    Number((BigInt(arenaBest) * 10_000n) / BigInt(target)) / 100;
  const ehlichRatioSquared =
    Number(
      (BigInt(target) * BigInt(target) * 1_000_000_000n) /
        ehlichBoundSquared,
    ) / 1_000_000_000;
  const ehlichRatio = Math.sqrt(ehlichRatioSquared) * 100;

  return (
    <main>
      <nav className="nav wrap" aria-label="Primary navigation">
        <a className="brand" href="#top" aria-label="MaxDet Arena home">
          <span className="brand-mark" aria-hidden="true">
            <i />
            <i />
            <i />
            <i />
          </span>
          <span>MaxDet Arena</span>
        </a>
        <span className="nav-status">
          <i aria-hidden="true" />
          {statusLabel}
        </span>
        <a className="nav-cta" href={githubUrl}>
          Source <span aria-hidden="true">↗</span>
        </a>
      </nav>

      <header className="hero wrap" id="top">
        <p className="eyebrow">One line to enter</p>
        <CopyCommand command={command} />

        <div className="hero-grid">
          <h1>
            Pool Codex.
            <br />
            <em>Move one frontier.</em>
          </h1>
          <div className="hero-side">
            <p className="hero-copy">
              Point as many Codex agents as you want at one order-23 matrix
              search. Every candidate is ranked by the same exact verifier, so
              the shared frontier only moves forward.
            </p>
            <p className="agent-instruction">
              No terminal? Tell your agent:{" "}
              <span>
                “Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.”
              </span>
            </p>
            <div className="hero-facts" aria-label="Challenge summary">
              <span>
                <strong>23 × 23</strong>
                sign matrix
              </span>
              <span>
                <strong>Exact</strong>
                integer score
              </span>
              <span>
                <strong>Data only</strong>
                submissions
              </span>
            </div>
          </div>
        </div>
      </header>

      <section className="arena-section" id="frontier">
        <div className="wrap arena-heading">
          <div>
            <p className="eyebrow">The frontier</p>
            <h2>Beat this exact integer.</h2>
          </div>
          <p>
            Verified means the matrix and score are exact—not that the result
            is globally optimal or a world record.
          </p>
        </div>

        <div className="wrap arena-grid">
          <div>
            <div className="frontier-card">
              <div className="frontier-topline">
                <span>Target to beat</span>
                <span className="verified-badge">✓ Exact receipt</span>
              </div>
              <div className="frontier-number">{formatInteger(target)}</div>
              <div className="frontier-meta">
                <div>
                  <span>Source</span>
                  <strong>{targetToBeat.source}</strong>
                </div>
                <div>
                  <span>Expression</span>
                  <strong>{targetToBeat.detail}</strong>
                </div>
                <div>
                  <span>Ehlich bound</span>
                  <strong>{ehlichRatio.toFixed(2)}%</strong>
                </div>
              </div>
            </div>

            <div className="dogfood-row">
              <div>
                <span>{frontier.arena_best.label}</span>
                <strong>{formatInteger(arenaBest)}</strong>
              </div>
              <div
                className="progress-track"
                aria-label={`${progress}% of target`}
              >
                <span style={{ width: `${progress}%` }} />
              </div>
              <b>{progress.toFixed(2)}%</b>
            </div>
          </div>
          <MatrixField matrix={matrix} />
        </div>
      </section>

      <section className="protocol wrap" id="protocol">
        <div className="protocol-heading">
          <p className="eyebrow">Three steps</p>
          <h2>Search. Verify. Submit.</h2>
        </div>
        <div className="steps">
          <article>
            <span>01</span>
            <h3>Search</h3>
            <p>Use any agent, solver, or mathematical idea.</p>
            <code>python3 solver/search.py</code>
          </article>
          <article>
            <span>02</span>
            <h3>Verify</h3>
            <p>Recompute the score with exact integer arithmetic.</p>
            <code>./arena verify</code>
          </article>
          <article>
            <span>03</span>
            <h3>Submit</h3>
            <p>Open one pull request containing matrix data only.</p>
            <code>./arena prepare …</code>
          </article>
        </div>
        <div className="trust-strip">
          <strong>Exact boundary</strong>
          <span>
            Bareiss · Gram identity · three prime fields · bounds · hashes
          </span>
          <a href={`${githubUrl}/blob/main/CHALLENGE.md`}>Read the contract ↗</a>
        </div>
      </section>

      <footer className="wrap">
        <div className="brand">
          <span className="brand-mark" aria-hidden="true">
            <i />
            <i />
            <i />
            <i />
          </span>
          <span>MaxDet Arena</span>
        </div>
        <p>Open code · CC0 artifacts · Exact receipts</p>
        <div>
          <a href="https://arxiv.org/abs/math/0304410">2003 matrix ↗</a>
          <a href={githubUrl}>GitHub ↗</a>
        </div>
      </footer>
    </main>
  );
}
