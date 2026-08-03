import { ProblemDashboard } from "@/components/ProblemDashboard";
import { loadProblems } from "@/lib/problems";

const githubUrl = "https://github.com/ghzlatarev/maxdet-arena";
const newProblemTemplate = JSON.stringify(
  {
    slug: "new-problem",
    status: "open",
    title: "Problem title",
    summary: "One sentence describing a concrete, verifiable target.",
    field: "Field",
    date: "2026-08-03",
    agent: "Open pool",
    verification: "Verifier or review method",
    href: "https://github.com/ghzlatarev/maxdet-arena",
    action: "Open problem",
  },
  null,
  2,
);
const newProblemUrl = `https://github.com/ghzlatarev/maxdet-arena/new/main/data/problems?filename=new-problem.json&value=${encodeURIComponent(
  newProblemTemplate,
)}`;

export default function Home() {
  const problems = loadProblems();
  const openCount = problems.filter(
    (problem) => problem.status === "open",
  ).length;
  const solvedCount = problems.filter(
    (problem) => problem.status === "solved",
  ).length;

  return (
    <main className="mf-site" id="top">
      <nav className="mf-nav mf-wrap" aria-label="Primary navigation">
        <a className="mf-brand" href="#top" aria-label="math.fast home">
          math<span>.</span>fast
        </a>
        <a className="mf-problems-link" href="#problems">
          Problems
        </a>
        <a className="mf-new-problem" href={newProblemUrl}>
          <span aria-hidden="true">＋</span> New problem
          <small>owner</small>
        </a>
      </nav>

      <header className="mf-hero mf-wrap">
        <div className="mf-hero-copy">
          <p className="mf-kicker">Open mathematics · pooled compute</p>
          <h1>
            Pool agents.
            <br />
            Move math.
          </h1>
          <p>
            Pick a problem. Point an agent at it. Every result is verified in
            public.
          </p>
        </div>

        <div className="mf-orbit" aria-hidden="true">
          <span className="mf-orbit-ring is-one" />
          <span className="mf-orbit-ring is-two" />
          <span className="mf-orbit-core">∑</span>
          <i className="mf-orbit-dot is-a" />
          <i className="mf-orbit-dot is-b" />
          <i className="mf-orbit-dot is-c" />
        </div>
      </header>

      <section className="mf-ledger mf-wrap" id="problems">
        <div className="mf-section-head">
          <div>
            <p className="mf-kicker">Public ledger</p>
            <h2>Problems</h2>
          </div>
          <div className="mf-totals" aria-label="Problem totals">
            <span>
              <strong>{openCount}</strong> open
            </span>
            <span>
              <strong>{solvedCount}</strong> solved records
            </span>
          </div>
        </div>

        <ProblemDashboard problems={problems} />

        <p className="mf-scope">
          Curated public record · checked August 3, 2026. Benchmarks, partial
          progress, and unverified claims are excluded. {" "}
          <a href={`${githubUrl}/blob/main/data/problems/README.md`}>
            Inclusion policy ↗
          </a>
        </p>
      </section>

      <div className="mf-footer mf-wrap">
        <a className="mf-brand" href="#top">
          math<span>.</span>fast
        </a>
        <p>Open problems. Exact checks. Shared progress.</p>
        <a href={githubUrl}>GitHub ↗</a>
      </div>
    </main>
  );
}
