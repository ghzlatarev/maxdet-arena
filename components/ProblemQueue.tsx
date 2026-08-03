"use client";

import Link from "next/link";
import { useState } from "react";
import { ProblemProposalTrigger } from "@/components/ProblemProposal";
import type { Problem } from "@/lib/problems";

type QueueView = "open" | "solved";

function ProblemAction({ href, label }: { href: string; label: string }) {
  const content = (
    <>
      {label} <span aria-hidden="true">→</span>
    </>
  );
  if (href.startsWith("/")) {
    return (
      <Link className="mf-open-action" href={href}>
        {content}
      </Link>
    );
  }
  return (
    <a className="mf-open-action" href={href}>
      {content}
    </a>
  );
}

export function ProblemQueue({ problems }: { problems: Problem[] }) {
  const [view, setView] = useState<QueueView>("open");
  const open = problems.filter((problem) => problem.status === "open");
  const solved = problems.filter((problem) => problem.status === "solved");

  return (
    <section className="mf-open-board" aria-labelledby="problem-queue-title">
      <h2 className="sr-only" id="problem-queue-title">
        Problem queue
      </h2>
      <header>
        <div className="mf-queue-tabs" aria-label="Problem views">
          <button
            aria-pressed={view === "open"}
            onClick={() => setView("open")}
            type="button"
          >
            Open problems
          </button>
          <button
            aria-pressed={view === "solved"}
            onClick={() => setView("solved")}
            type="button"
          >
            Solved
          </button>
        </div>
        <span aria-live="polite">
          {view === "open" ? `${open.length} active` : `${solved.length} records`}
        </span>
      </header>

      {view === "open" ? (
        <div aria-label="Open problem list" className="mf-open-list" role="list">
          {open.map((problem, index) => (
            <article key={problem.slug} role="listitem">
              <span className="mf-open-number" aria-hidden="true">
                {String(index + 1).padStart(2, "0")}
              </span>
              <div className="mf-open-summary">
                <div className="mf-open-topline">
                  <span>
                    <i aria-hidden="true" /> Open
                  </span>
                  <span>{problem.field}</span>
                </div>
                <h3>{problem.title}</h3>
                <p>{problem.summary}</p>
              </div>
              <dl className="mf-open-meta">
                <div>
                  <dt>Verification</dt>
                  <dd>{problem.verification}</dd>
                </div>
                <div>
                  <dt>Pool</dt>
                  <dd>{problem.agent}</dd>
                </div>
              </dl>
              <ProblemAction href={problem.href} label={problem.action} />
            </article>
          ))}
          <ProblemProposalTrigger position={open.length + 1} />
        </div>
      ) : (
        <div
          aria-label="Solved problem records"
          className="mf-solved-queue"
          role="list"
        >
          {solved.map((problem, index) => (
            <article key={problem.slug} role="listitem">
              <span className="mf-solved-number" aria-hidden="true">
                {String(index + 1).padStart(2, "0")}
              </span>
              <div className="mf-solved-summary">
                <div>
                  <span>Solved</span>
                  <span>{problem.field}</span>
                </div>
                <h3>{problem.title}</h3>
              </div>
              <span className="mf-solved-check">{problem.verification}</span>
              <ProblemAction href={problem.href} label="View record" />
            </article>
          ))}
        </div>
      )}
    </section>
  );
}
