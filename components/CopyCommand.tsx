"use client";

import { useEffect, useState } from "react";

type CopyCommandProps = {
  command: string;
};

export function CopyCommand({ command }: CopyCommandProps) {
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    if (!copied) return;
    const timeout = window.setTimeout(() => setCopied(false), 1800);
    return () => window.clearTimeout(timeout);
  }, [copied]);

  async function copy() {
    try {
      await navigator.clipboard.writeText(command);
    } catch {
      const fallback = document.createElement("textarea");
      fallback.value = command;
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
    <div className="command-shell">
      <div className="command-label">
        <span className="terminal-dot" aria-hidden="true" />
        One line to enter the arena
      </div>
      <div className="command-row">
        <code>{command}</code>
        <button
          type="button"
          onClick={copy}
          aria-label="Copy terminal command"
        >
          {copied ? "Copied" : "Copy"}
        </button>
        <span className="sr-only" role="status" aria-live="polite">
          {copied ? "Terminal command copied." : ""}
        </span>
      </div>
    </div>
  );
}
