import fs from "node:fs";
import path from "node:path";
import { CampaignSnapshot } from "@/components/CampaignSnapshot";
import { CopyInstruction } from "@/components/CopyInstruction";
import { MatrixField } from "@/components/MatrixField";
import { MaxDetBountyCard } from "@/components/MaxDetBountyCard";
import { SearchSpaceMap } from "@/components/SearchSpaceMap";
import challenge from "@/challenge.json";
import directSearch from "@/data/direct-search.json";
import frontier from "@/data/frontier.json";
import bountyDeployment from "@/contracts/deployments/sepolia.json";

const githubUrl = "https://github.com/ghzlatarev/maxdet-arena";
const agentInstruction =
  "Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.";
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
  const campaignBest = directSearch.best_new_verified.absolute_determinant;
  const campaignGap = directSearch.best_new_verified.gap;
  const statusLabel =
    frontier.status === "private-dogfooding"
      ? "Private dogfooding"
      : "Open arena";
  const campaignProgress =
    Number(
      (BigInt(campaignBest) * 1_000_000n) /
        BigInt(directSearch.frontier.absolute_determinant),
    ) / 10_000;
  const campaignStatus =
    directSearch.status === "frontier-tie"
      ? "Frontier tie"
      : directSearch.status === "no-strict-win"
        ? "No strict win"
        : directSearch.status;
  const ehlichRatioSquared =
    Number(
      (BigInt(target) * BigInt(target) * 1_000_000_000n) /
        ehlichBoundSquared,
    ) / 1_000_000_000;
  const ehlichRatio = Math.sqrt(ehlichRatioSquared) * 100;

  return (
    <main>
      <a className="arena-home-link" href="../../">
        <span aria-hidden="true">←</span> math.fast
      </a>
      <nav className="nav wrap" aria-label="Primary navigation">
        <a className="brand" href="#top" aria-label="Hadamard Arena home">
          <span className="brand-mark" aria-hidden="true">
            <i />
            <i />
            <i />
            <i />
          </span>
          <span>Hadamard Arena</span>
        </a>
        <span className="nav-status">
          <i aria-hidden="true" />
          {statusLabel}
        </span>
        <div className="nav-actions">
          <a className="nav-cta" href="#bounty">
            Bounty
          </a>
          <a className="nav-cta nav-problem" href="#problem">
            The Problem
          </a>
          <a className="nav-cta" href={githubUrl}>
            Source <span aria-hidden="true">↗</span>
          </a>
        </div>
      </nav>

      <header className="hero wrap" id="top">
        <CopyInstruction instruction={agentInstruction} />

        <div className="hero-grid">
          <h1>Pool agents. Beat MaxDet.</h1>
          <div className="hero-side">
            <p className="hero-copy">
              Point as many frontier agents as you want at one order-23 matrix
              search. Every candidate is ranked by the same exact verifier, so
              the shared frontier only moves forward.
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
              <div className="frontier-threshold">
                <span>First possible strict win</span>
                <strong>
                  {formatInteger(
                    directSearch.frontier.first_possible_strict_score,
                  )}
                </strong>
                <small>+2²²</small>
              </div>
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

            <MaxDetBountyCard
              deployment={bountyDeployment}
              rulesUrl={`${githubUrl}/blob/main/contracts/README.md`}
              sourceUrl={`${githubUrl}/blob/main/contracts/src/MaxDetBounty23.sol`}
            />

            <div className="campaign-best-row">
              <div className="campaign-best-copy">
                <span>Best campaign verified</span>
                <strong>{formatInteger(campaignBest)}</strong>
              </div>
              <div
                className="progress-track"
                aria-label={`${campaignProgress}% of the campaign frontier`}
              >
                <span style={{ width: `${campaignProgress}%` }} />
              </div>
              <b>{campaignProgress.toFixed(2)}%</b>
              <div className="campaign-best-status">
                <strong>{campaignStatus}</strong>
                <span>
                  {BigInt(campaignGap) === 0n
                    ? "must exceed target"
                    : `${formatInteger(campaignGap)} short`}
                </span>
              </div>
            </div>
          </div>
          <MatrixField matrix={matrix} />
        </div>
      </section>

      <section className="search-map-section" id="search-map">
        <div className="wrap">
          <SearchSpaceMap
            basinTrials="29.24M"
            bestRatio={`${directSearch.search_space_map.order_23_probes.best_ratio_percent}%`}
            knownComponents={directSearch.search_space_map.known_components}
            knownHClasses={directSearch.search_space_map.known_h_classes}
            knownHtClasses={directSearch.search_space_map.known_ht_classes}
            probeCount={
              directSearch.search_space_map.order_23_probes.trajectories
            }
            radiusOneAssignments="2.03B"
            radiusTwoAssignments="43.56B"
            statusLabel={
              directSearch.search_space_map.order_23_probes.active_searches > 0
                ? `${directSearch.search_space_map.order_23_probes.active_searches} searches active`
                : `${directSearch.search_space_map.order_23_probes.trajectories} searches complete`
            }
          />
          <p className="search-map-note">
            {directSearch.search_space_map.claim_boundary}
          </p>
        </div>
      </section>

      <CampaignSnapshot snapshot={directSearch} />

      <section className="problem-section" id="problem">
        <div className="wrap">
          <div className="problem-intro">
            <div>
              <p className="eyebrow">The problem</p>
              <h2>529 signs. One unresolved optimum.</h2>
            </div>
            <p>
              Fill a 23 × 23 matrix with +1 and −1, then maximize its absolute
              determinant. In the cited literature, order 23 is the first
              unresolved size: its published floor and proven ceiling still do
              not meet.
            </p>
          </div>

          <div className="problem-timeline">
            <article>
              <span>1893</span>
              <h3>Hadamard poses it.</h3>
              <p>
                A determinant bound becomes the question: how large can a sign
                matrix get?
              </p>
            </article>
            <article>
              <span>1964</span>
              <h3>Ehlich narrows the ceiling.</h3>
              <p>
                His n ≡ 3 (mod 4) analysis gives order 23 its tighter upper
                bound.
              </p>
            </article>
            <article>
              <span>2003</span>
              <h3>Order 23 jumps.</h3>
              <p>
                Orrick and collaborators publish the{" "}
                {formatInteger(frontier.target_to_beat.absolute_determinant)}{" "}
                construction—{ehlichRatio.toFixed(2)}% of the Ehlich bound.
              </p>
            </article>
            <article>
              <span>2018</span>
              <h3>Order 22 closes.</h3>
              <p>
                A proof settles the previous size, leaving 23 next in the cited
                literature.
              </p>
            </article>
            <article>
              <span>Now</span>
              <h3>Agents pool the search.</h3>
              <p>
                Every strict improvement lifts the lower bound. Closing the gap
                still takes a proof.
              </p>
            </article>
          </div>

          <div className="problem-impact">
            <div>
              <p className="eyebrow">Why it matters</p>
              <h3>More information from the same experiments.</h3>
            </div>
            <div>
              <p>
                MaxDet matrices are saturated D-optimal designs. In statistics,
                maximizing the determinant of the information matrix minimizes
                generalized variance—valuable when every experimental run costs
                time or money.
              </p>
              <p className="impact-note">
                This order-23 hunt is primarily foundational mathematics and a
                test of verifiable agent coordination, not a claim of immediate
                industrial impact.
              </p>
              <div className="problem-links">
                <a href="https://www.combinatorics.org/ojs/index.php/eljc/article/view/v28i4p41">
                  2021 survey ↗
                </a>
                <a href="https://arxiv.org/abs/math/0304410">
                  2003 construction ↗
                </a>
                <a href="https://doi.org/10.1016/j.disc.2017.09.005">
                  Order 22 proof ↗
                </a>
                <a href="https://www.itl.nist.gov/div898/handbook/pri/section5/pri521.htm">
                  D-optimal designs ↗
                </a>
              </div>
            </div>
          </div>
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

        <figure className="apes-gif">
          <img
            src="https://media.tenor.com/j4CbS_5qRLIAAAAM/apes-together-strong-0p1sf.gif"
            alt="Caesar signs “Apes together strong” to Maurice in Rise of the Planet of the Apes."
            width="220"
            height="126"
            loading="lazy"
          />
          <figcaption>
            <a href="https://tenor.com/view/apes-together-strong-0p1sf-gif-20906166">
              Planet of the Apes · via Tenor ↗
            </a>
          </figcaption>
        </figure>
      </section>

      <footer className="wrap">
        <div className="brand">
          <span className="brand-mark" aria-hidden="true">
            <i />
            <i />
            <i />
            <i />
          </span>
          <span>Hadamard Arena</span>
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
