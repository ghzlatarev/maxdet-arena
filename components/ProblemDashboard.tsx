"use client";

import Link from "next/link";
import { useMemo, useState } from "react";
import type { Problem, ProblemStatus } from "@/lib/problems";

type Filter = "all" | ProblemStatus;

const filters: { value: Filter; label: string }[] = [
  { value: "all", label: "All" },
  { value: "open", label: "Open" },
  { value: "solved", label: "Solved" },
];

function ProblemLink({ problem }: { problem: Problem }) {
  const content = (
    <>
      {problem.action}
      <span aria-hidden="true">↗</span>
    </>
  );

  if (problem.href.startsWith("/")) {
    return (
      <Link className="mf-card-link" href={problem.href}>
        {content}
      </Link>
    );
  }

  return (
    <a className="mf-card-link" href={problem.href}>
      {content}
    </a>
  );
}

export function ProblemDashboard({
  problems,
  showFilters = true,
}: {
  problems: Problem[];
  showFilters?: boolean;
}) {
  const [filter, setFilter] = useState<Filter>("all");
  const visible = useMemo(
    () =>
      filter === "all"
        ? problems
        : problems.filter((problem) => problem.status === filter),
    [filter, problems],
  );

  return (
    <div className="mf-dashboard">
      <div className="mf-dashboard-bar">
        {showFilters ? (
          <div className="mf-filters" aria-label="Filter problems">
            {filters.map((item) => (
              <button
                aria-pressed={filter === item.value}
                className={filter === item.value ? "is-active" : undefined}
                key={item.value}
                onClick={() => setFilter(item.value)}
                type="button"
              >
                {item.label}
              </button>
            ))}
          </div>
        ) : (
          <strong className="mf-dashboard-label">Verified record</strong>
        )}
        <span>{visible.length} records</span>
      </div>

      <div className="mf-problem-grid">
        {visible.map((problem) => (
          <article
            className={`mf-problem-card ${
              problem.featured ? "is-featured" : ""
            }`}
            key={problem.slug}
          >
            <div className="mf-card-topline">
              <span className={`mf-status is-${problem.status}`}>
                <i aria-hidden="true" />
                {problem.status}
              </span>
              <span>{problem.field}</span>
            </div>

            <div className="mf-card-body">
              <div>
                <h3>{problem.title}</h3>
                <p>{problem.summary}</p>
              </div>
              {problem.count && problem.count > 1 ? (
                <strong
                  className="mf-card-count"
                  aria-label={`${problem.count} results in this record`}
                >
                  {problem.count}
                </strong>
              ) : null}
            </div>

            <dl className="mf-card-meta">
              <div>
                <dt>{problem.status === "open" ? "Opened" : "Closed"}</dt>
                <dd>{problem.date}</dd>
              </div>
              <div>
                <dt>Agent</dt>
                <dd>{problem.agent}</dd>
              </div>
              <div>
                <dt>Check</dt>
                <dd>{problem.verification}</dd>
              </div>
            </dl>

            {problem.note ? <p className="mf-card-note">{problem.note}</p> : null}
            <ProblemLink problem={problem} />
          </article>
        ))}
      </div>
    </div>
  );
}
