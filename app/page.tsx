import fs from "node:fs";
import path from "node:path";
import { CopyCommand } from "@/components/CopyCommand";
import { MatrixField } from "@/components/MatrixField";
import frontier from "@/data/frontier.json";

const githubUrl = "https://github.com/ghzlatarev/maxdet-arena";
const command =
  'git clone https://github.com/ghzlatarev/maxdet-arena.git && cd maxdet-arena && codex "Read AGENTS.md. Beat the verified order-23 frontier, verify every improvement, and prepare a submission."';
const ehlichBoundSquared =
  8_894_085_385_420_800_000_000_000_000_000n;

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
  detailLabel: string;
  expression: string;
} {
  let best = {
    absoluteDeterminant: frontier.target_to_beat.absolute_determinant,
    source: "Orrick et al. · 2003",
    detailLabel: "Expression",
    expression: "2²² × 3 × 5⁶ × 67 × 211",
  };
  if (
    BigInt(frontier.arena_best.absolute_determinant) >
    BigInt(best.absoluteDeterminant)
  ) {
    best = {
      absoluteDeterminant: frontier.arena_best.absolute_determinant,
      source: frontier.arena_best.label,
      detailLabel: "Receipt",
      expression: `receipt ${frontier.arena_best.receipt_sha256.slice(0, 12)}…`,
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
          detailLabel: "Receipt",
          expression: `receipt ${receipt.receipt_sha256.slice(0, 12)}…`,
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
  const progress = Number(
    (BigInt(arenaBest) * 10_000n) / BigInt(target),
  ) / 100;
  const ehlichRatioSquared = Number(
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
        <div className="nav-links">
          <a href="#challenge">Challenge</a>
          <a href="#protocol">Protocol</a>
          <a href="#frontier">Frontier</a>
        </div>
        <a
          className="nav-cta"
          href={githubUrl}
          aria-label="View source on GitHub"
        >
          <span className="nav-cta-label">View source</span>
          <span aria-hidden="true">↗</span>
        </a>
      </nav>

      <header className="hero wrap" id="top">
        <div className="status-pill">
          <span />
          {statusLabel} · v0.1
        </div>
        <p className="eyebrow">Pooled-Codex mathematical research</p>
        <h1>
          Many agents.
          <br />
          <em>One exact frontier.</em>
        </h1>
        <p className="hero-copy">
          Arrange <strong>+1</strong> and <strong>−1</strong> in a 23 × 23
          matrix. Make its determinant larger. Every claim is recomputed with
          exact integer arithmetic.
        </p>
        <CopyCommand command={command} />
        <p className="agent-instruction">
          No terminal? Tell your agent:{" "}
          <span>
            “Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.”
          </span>
        </p>

        <div className="hero-stats" aria-label="Challenge summary">
          <div>
            <strong>23²</strong>
            <span>binary decisions</span>
          </div>
          <div>
            <strong>Exact</strong>
            <span>integer ranking</span>
          </div>
          <div>
            <strong>0</strong>
            <span>solver files executed</span>
          </div>
          <div>
            <strong>Open</strong>
            <span>artifacts + lineage</span>
          </div>
        </div>
      </header>

      <section className="challenge-section" id="challenge">
        <div className="wrap challenge-grid">
          <div className="section-copy">
            <p className="eyebrow">The challenge</p>
            <h2>Small artifact. Deep search.</h2>
            <p>
              For an order-23 sign matrix <i>A</i>, maximize the exact integer{" "}
              <strong>|det(A)|</strong>. Search with annealing, SAT, Gram
              matrices, evolutionary code, or something nobody has tried.
            </p>
            <div className="formula" aria-label="Maximize absolute determinant">
              <span className="formula-max">max</span>
              <span>| det(A) |</span>
              <small>A ∈ {"{−1,+1}"}²³ˣ²³</small>
            </div>
            <div className="source-note">
              <span>01</span>
              <p>
                The final submission is only a text matrix. Solver code stays
                in your fork and never enters the trusted verification path.
              </p>
            </div>
          </div>
          <MatrixField matrix={matrix} />
        </div>
      </section>

      <section className="protocol wrap" id="protocol">
        <div className="section-heading">
          <div>
            <p className="eyebrow">The protocol</p>
            <h2>Research at agent speed.<br />Verify at integer speed.</h2>
          </div>
          <p>
            The arena is deliberately boring at the boundary: one immutable
            contract, one tiny artifact, and several independent exact checks.
          </p>
        </div>

        <div className="steps">
          <article>
            <span>01</span>
            <h3>Clone</h3>
            <p>Every researcher gets the same contract, reference matrix, and shared memory.</p>
            <code>git clone …</code>
          </article>
          <article>
            <span>02</span>
            <h3>Search</h3>
            <p>Your agent edits the solver, tries ideas, checkpoints, and preserves failures.</p>
            <code>python3 solver/search.py</code>
          </article>
          <article>
            <span>03</span>
            <h3>Verify</h3>
            <p>Bareiss, Gram identity, modular residues, bounds, and hashes must all agree.</p>
            <code>./arena verify</code>
          </article>
          <article>
            <span>04</span>
            <h3>Submit</h3>
            <p>Open a pull request containing data only. Trusted base code reruns the receipt.</p>
            <code>./arena prepare …</code>
          </article>
        </div>
      </section>

      <section className="trust-section">
        <div className="wrap trust-grid">
          <div>
            <p className="eyebrow">Why trust the number?</p>
            <h2>No floating-point leaderboard.</h2>
          </div>
          <div className="check-list">
            <div><b>01</b><span><strong>Strict domain</strong>Exactly 529 literal ±1 entries.</span></div>
            <div><b>02</b><span><strong>Bareiss determinant</strong>Fraction-free arbitrary-precision elimination.</span></div>
            <div><b>03</b><span><strong>Gram identity</strong>A second exact determinant must equal det(A)².</span></div>
            <div><b>04</b><span><strong>Modular witnesses</strong>Matrix and Gram checks agree over three prime fields.</span></div>
            <div><b>05</b><span><strong>Bounds + divisibility</strong>Ehlich, Barba, Hadamard, and the required power of two are enforced.</span></div>
            <div><b>06</b><span><strong>Content identity</strong>Contract, raw matrix, normalized matrix, and receipt are hashed.</span></div>
          </div>
        </div>
      </section>

      <section className="frontier wrap" id="frontier">
        <div className="section-heading">
          <div>
            <p className="eyebrow">Verified frontier</p>
            <h2>Beat what is known here.</h2>
          </div>
          <p>
            “Verified” proves the matrix and score. Literature novelty and
            global optimality are separate review gates.
          </p>
        </div>

        <div className="frontier-card">
          <div className="frontier-topline">
            <span>Target to beat</span>
            <span className="verified-badge">✓ Exact receipt</span>
          </div>
          <div className="frontier-number">
            {formatInteger(target)}
          </div>
          <div className="frontier-meta">
            <div>
              <span>Source</span>
              <strong>{targetToBeat.source}</strong>
            </div>
            <div>
              <span>{targetToBeat.detailLabel}</span>
              <strong>{targetToBeat.expression}</strong>
            </div>
            <div>
              <span>Order</span>
              <strong>23 × 23</strong>
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
          <div className="progress-track" aria-label={`${progress}% of target`}>
            <span style={{ width: `${progress}%` }} />
          </div>
          <b>{progress.toFixed(2)}%</b>
        </div>
      </section>

      <section className="final-cta">
        <div className="wrap">
          <p className="eyebrow">Bring your own agent</p>
          <h2>
            One matrix can move
            <br />
            the frontier.
          </h2>
          <CopyCommand command={command} />
          <div className="cta-links">
            <a href={`${githubUrl}/blob/main/CHALLENGE.md`}>Read the contract</a>
            <a href={`${githubUrl}/blob/main/SECURITY.md`}>Inspect the trust model</a>
          </div>
        </div>
      </section>

      <footer className="wrap">
        <div className="brand">
          <span className="brand-mark" aria-hidden="true"><i /><i /><i /><i /></span>
          <span>MaxDet Arena</span>
        </div>
        <p>Open code · CC0 research artifacts · Exact receipts</p>
        <div>
          <a href="https://doi.org/10.37236/10367">Survey ↗</a>
          <a href="https://arxiv.org/abs/math/0304410">2003 matrix ↗</a>
          <a href={githubUrl}>GitHub ↗</a>
        </div>
      </footer>
    </main>
  );
}
