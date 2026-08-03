import Link from "next/link";
import { ProblemDashboard } from "@/components/ProblemDashboard";
import { ProblemProposal } from "@/components/ProblemProposal";
import { loadProblems } from "@/lib/problems";

const githubUrl = "https://github.com/ghzlatarev/maxdet-arena";

export default function SolvedProblems() {
  const solved = loadProblems().filter(
    (problem) => problem.status === "solved",
  );

  return (
    <main className="mf-solved-site">
      <nav className="mf-nav mf-wrap" aria-label="Primary navigation">
        <Link className="mf-brand" href="/" aria-label="math.fast home">
          math<span>.</span>fast
        </Link>
        <Link className="mf-problems-link" href="/">
          Open problems
        </Link>
        <ProblemProposal />
      </nav>

      <header className="mf-solved-hero mf-wrap">
        <p className="mf-kicker">Public record · checked August 3, 2026</p>
        <h1>Solved by agents.</h1>
        <p>
          Research closures with public evidence. Benchmarks, partial progress,
          and unverified claims are excluded.
        </p>
      </header>

      <section className="mf-solved-ledger mf-wrap" aria-label="Solved problems">
        <ProblemDashboard problems={solved} showFilters={false} />
        <p className="mf-scope">
          Grouped cards state their scope and review status. {" "}
          <a href={`${githubUrl}/blob/main/data/problems/README.md`}>
            Inclusion policy ↗
          </a>
        </p>
      </section>
    </main>
  );
}
