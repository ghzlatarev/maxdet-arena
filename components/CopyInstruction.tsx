"use client";

import { useEffect, useState } from "react";

type CopyInstructionProps = {
  instruction: string;
};

export function CopyInstruction({ instruction }: CopyInstructionProps) {
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    if (!copied) return;
    const timeout = window.setTimeout(() => setCopied(false), 1800);
    return () => window.clearTimeout(timeout);
  }, [copied]);

  async function copy() {
    try {
      await navigator.clipboard.writeText(instruction);
    } catch {
      const fallback = document.createElement("textarea");
      fallback.value = instruction;
      fallback.style.position = "fixed";
      fallback.style.opacity = "0";
      document.body.appendChild(fallback);
      fallback.select();
      const copiedWithFallback = document.execCommand("copy");
      fallback.remove();
      if (!copiedWithFallback) return;
    }
    setCopied(true);
  }

  return (
    <div className="agent-ask">
      <div className="agent-ask-label">Simply ask your agent:</div>
      <div className="agent-ask-row">
        <p>{instruction}</p>
        <button
          type="button"
          onClick={copy}
          aria-label="Copy agent instruction"
        >
          {copied ? "Copied" : "Copy"}
        </button>
        <span className="sr-only" role="status" aria-live="polite">
          {copied ? "Agent instruction copied." : ""}
        </span>
      </div>
    </div>
  );
}
