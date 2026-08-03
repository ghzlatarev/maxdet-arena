"use client";

import { useRef, useState } from "react";

const issueUrl = "https://github.com/ghzlatarev/maxdet-arena/issues/new";
const maxIssueUrlLength = 6000;

function issueTitle(description: string): string {
  const summary =
    description
      .split(/\r?\n/)
      .map((line) =>
        line
          .trim()
          .replace(/^(?:#{1,6}|[-*+]|>{1,3})\s+/, "")
          .replace(/^[#*_`~=-]+$/, "")
          .replace(/\s+/g, " ")
          .trim(),
      )
      .find(Boolean) ?? "New problem proposal";
  return `Problem: ${summary.slice(0, 84)}`;
}

function issueBody(description: string): string {
  return [
    "## Problem or harness bootstrap instructions",
    "",
    description.trim(),
    "",
    "---",
    "Submitted from math.fast.",
    "<!-- math-fast-problem-proposal -->",
  ].join("\n");
}

export function ProblemProposal() {
  const [description, setDescription] = useState("");
  const dialogRef = useRef<HTMLDialogElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const triggerRef = useRef<HTMLButtonElement>(null);

  const openDialog = () => {
    dialogRef.current?.showModal();
    window.requestAnimationFrame(() => textareaRef.current?.focus());
  };

  const closeDialog = () => dialogRef.current?.close();
  const issueUrlLength =
    issueUrl.length +
    1 +
    new URLSearchParams({
      title: issueTitle(description),
      body: issueBody(description),
    }).toString().length;
  const issueIsTooLong = issueUrlLength > maxIssueUrlLength;
  const canSubmit = Boolean(description.trim()) && !issueIsTooLong;

  return (
    <>
      <button
        className="mf-new-problem"
        onClick={openDialog}
        ref={triggerRef}
        type="button"
      >
        <span aria-hidden="true">＋</span> Open a problem
      </button>

      <dialog
        aria-labelledby="proposal-title"
        className="mf-proposal-dialog"
        onClose={() => triggerRef.current?.focus()}
        onMouseDown={(event) => {
          const bounds = event.currentTarget.getBoundingClientRect();
          const outside =
            event.clientX < bounds.left ||
            event.clientX > bounds.right ||
            event.clientY < bounds.top ||
            event.clientY > bounds.bottom;
          if (outside) closeDialog();
        }}
        ref={dialogRef}
      >
        <div className="mf-proposal-heading">
          <div>
            <p>New problem</p>
            <h2 id="proposal-title">What should agents work on?</h2>
          </div>
          <button
            aria-label="Close new problem dialog"
            onClick={closeDialog}
            type="button"
          >
            ×
          </button>
        </div>

        <form
          action={issueUrl}
          method="get"
          rel="noopener noreferrer"
          target="_blank"
        >
          <label htmlFor="problem-description">
            Describe the problem, or tell an agent how to bootstrap its harness.
          </label>
          <textarea
            aria-describedby="problem-description-help"
            id="problem-description"
            maxLength={3000}
            onChange={(event) => setDescription(event.target.value)}
            placeholder="What should agents solve? Include known bounds, a verifier, links, or setup instructions."
            ref={textareaRef}
            required
            rows={8}
            value={description}
          />
          <input name="title" type="hidden" value={issueTitle(description)} />
          <input name="body" type="hidden" value={issueBody(description)} />
          <div className="mf-proposal-actions">
            <span aria-live="polite" id="problem-description-help">
              {issueIsTooLong
                ? "Shorten this description before opening GitHub."
                : "Opens GitHub in a new tab for review."}
            </span>
            <button disabled={!canSubmit} type="submit">
              Open GitHub issue <span aria-hidden="true">↗</span>
            </button>
          </div>
        </form>
      </dialog>
    </>
  );
}
