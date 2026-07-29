"use client";

import { useEffect, useMemo, useState } from "react";

import styles from "./SearchSpaceMap.module.css";

export type SearchSpaceMapProps = {
  className?: string;
  basinTrials?: string;
  bestRatio?: string;
  knownComponents?: number;
  knownHClasses?: number;
  knownHtClasses?: number;
  probeCount?: number;
  radiusOneAssignments?: string;
  radiusTwoAssignments?: string;
  statusLabel?: string;
};

type LiveSearchArm = {
  id: string;
  label: string;
  engine: string;
  status: "active" | "complete";
  elapsed_seconds: number;
  budget_seconds: number;
  progress_percent: number;
  best_absolute_determinant: string;
  best_ratio_percent: number;
};

const LIVE_TELEMETRY_MAX_AGE_MS = 15_000;

type LiveSearchProgress = {
  updated_at: string;
  frontier_absolute_determinant: string;
  arms: LiveSearchArm[];
};

const historicalRoutePaths = [
  "M282 136 C434 85 632 79 807 119",
  "M270 326 C436 356 627 321 804 275",
  "M506 167 C625 184 721 193 812 197",
  "M632 321 C703 330 754 316 805 300",
] as const;

const routePaths = [
  "M676 115 C731 113 771 116 807 119",
  "M620 144 C697 192 754 240 804 275",
  "M681 146 C731 153 773 174 812 197",
  "M651 161 C706 205 760 259 805 300",
] as const;

function formatDuration(seconds: number) {
  const wholeSeconds = Math.max(0, Math.round(seconds));
  const minutes = Math.floor(wholeSeconds / 60);
  const remainder = wholeSeconds % 60;
  return `${minutes}:${remainder.toString().padStart(2, "0")}`;
}

function isLiveSearchProgress(value: unknown): value is LiveSearchProgress {
  if (!value || typeof value !== "object") {
    return false;
  }

  const candidate = value as Partial<LiveSearchProgress>;
  return (
    typeof candidate.updated_at === "string" &&
    typeof candidate.frontier_absolute_determinant === "string" &&
    Array.isArray(candidate.arms) &&
    candidate.arms.every(
      (arm) =>
        arm &&
        typeof arm.id === "string" &&
        typeof arm.label === "string" &&
        typeof arm.engine === "string" &&
        (arm.status === "active" || arm.status === "complete") &&
        typeof arm.elapsed_seconds === "number" &&
        typeof arm.budget_seconds === "number" &&
        typeof arm.progress_percent === "number" &&
        typeof arm.best_absolute_determinant === "string" &&
        typeof arm.best_ratio_percent === "number",
    )
  );
}

const alphaNodes = [
  [174, 115],
  [204, 105],
  [235, 113],
  [258, 136],
  [249, 165],
  [219, 178],
  [187, 169],
  [164, 146],
  [211, 139],
  [229, 147],
] as const;

const betaNodes = [
  [172, 316],
  [206, 306],
  [219, 338],
  [184, 348],
] as const;

const gammaNodes = [
  [421, 154],
  [455, 145],
  [469, 178],
  [433, 190],
] as const;

const deltaNodes = [
  [529, 293],
  [563, 284],
  [594, 300],
  [603, 332],
  [570, 349],
  [537, 333],
] as const;

const epsilonNodes = [
  [616, 111],
  [647, 98],
  [676, 115],
  [681, 146],
  [651, 161],
  [620, 144],
] as const;

function BasinNodes({
  nodes,
  accent,
}: {
  nodes: readonly (readonly [number, number])[];
  accent: "cyan" | "violet" | "amber" | "rose" | "blue";
}) {
  return (
    <g className={styles.nodes} data-accent={accent} aria-hidden="true">
      {nodes.map(([x, y], index) => (
        <g key={`${x}-${y}`}>
          <circle className={styles.nodeHalo} cx={x} cy={y} r="8" />
          <circle className={styles.node} cx={x} cy={y} r="3.2" />
          {index === 0 ? (
            <circle className={styles.seedNode} cx={x} cy={y} r="6.4" />
          ) : null}
        </g>
      ))}
    </g>
  );
}

export function SearchSpaceMap({
  className,
  basinTrials = "29.24M",
  bestRatio = "96.96%",
  knownComponents = 5,
  knownHClasses = 30,
  knownHtClasses = 26,
  probeCount = 4,
  radiusOneAssignments = "2.03B",
  radiusTwoAssignments = "43.56B",
  statusLabel = "2 exact refinements active",
}: SearchSpaceMapProps) {
  const rootClassName = [styles.root, className].filter(Boolean).join(" ");
  const [liveProgress, setLiveProgress] =
    useState<LiveSearchProgress | null>(null);
  const [telemetryClock, setTelemetryClock] = useState(0);

  useEffect(() => {
    let mounted = true;
    const progressUrl = new URL("./search-progress.json", window.location.href);

    const refresh = async () => {
      try {
        const response = await fetch(`${progressUrl}?t=${Date.now()}`, {
          cache: "no-store",
        });
        if (!response.ok) {
          return;
        }
        const payload: unknown = await response.json();
        if (mounted && isLiveSearchProgress(payload)) {
          setLiveProgress(payload);
          setTelemetryClock(Date.now());
        }
      } catch {
        // The checked-in campaign snapshot remains visible when live telemetry
        // is unavailable (for example, on an offline static export).
      }
    };

    void refresh();
    const interval = window.setInterval(() => void refresh(), 2_000);
    return () => {
      mounted = false;
      window.clearInterval(interval);
    };
  }, []);

  const liveArms = liveProgress?.arms.slice(0, routePaths.length) ?? [];
  const telemetryAge = liveProgress
    ? telemetryClock - Date.parse(liveProgress.updated_at)
    : Number.POSITIVE_INFINITY;
  const telemetryIsFresh =
    Number.isFinite(telemetryAge) &&
    telemetryAge >= 0 &&
    telemetryAge <= LIVE_TELEMETRY_MAX_AGE_MS;
  const liveBestRatio = useMemo(() => {
    if (liveArms.length === 0) {
      return bestRatio;
    }
    return `${Math.max(...liveArms.map((arm) => arm.best_ratio_percent)).toFixed(2)}%`;
  }, [bestRatio, liveArms]);
  const activeSearches = telemetryIsFresh
    ? liveArms.filter((arm) => arm.status === "active").length
    : 0;
  const resolvedStatus =
    liveArms.length > 0
      ? !telemetryIsFresh
        ? `${liveArms.length}-probe telemetry snapshot`
        : activeSearches > 0
        ? `${activeSearches} searches active`
        : `${liveArms.length} searches complete`
      : statusLabel;

  return (
    <section
      className={rootClassName}
      aria-label="Exact search-space map for the order-22 and order-23 MaxDet campaign"
    >
      <header className={styles.header}>
        <div>
          <p className={styles.eyebrow}>Search-space projection · 2^529</p>
          <h3>Where the search is looking.</h3>
        </div>
        <dl className={styles.metrics} aria-label="Mapped search-space totals">
          <div>
            <dt>Basins</dt>
            <dd>{knownComponents}</dd>
          </div>
          <div>
            <dt>H classes</dt>
            <dd>{knownHClasses}</dd>
          </div>
          <div>
            <dt>Best probe</dt>
            <dd>{liveBestRatio}</dd>
          </div>
        </dl>
      </header>

      <div className={styles.mapFrame}>
        <svg
          className={styles.map}
          viewBox="0 0 920 460"
          role="img"
          aria-label="Five known order-22 maximal determinant basins with exact radius-one and radius-two coverage, and animated order-23 probes approaching the current frontier"
          preserveAspectRatio="xMidYMid meet"
        >
          <title>MaxDet search-space atlas</title>
          <desc>
            Five disconnected order-22 components contain 30 H-equivalence
            classes and 26 HT-equivalence classes. Concentric rings around the
            two published seeds show exact radius-one and radius-two coverage.
            Three additional basins were recovered by randomized harvesting.
            Deterministic moving probes travel toward the order-23 frontier.
          </desc>

          <g className={styles.grid} aria-hidden="true">
            {Array.from({ length: 17 }, (_, index) => (
              <line
                key={`vertical-${index}`}
                x1={32 + index * 52}
                x2={32 + index * 52}
                y1="30"
                y2="430"
              />
            ))}
            {Array.from({ length: 8 }, (_, index) => (
              <line
                key={`horizontal-${index}`}
                x1="32"
                x2="888"
                y1={45 + index * 52}
                y2={45 + index * 52}
              />
            ))}
          </g>

          <path
            className={styles.unmappedField}
            d="M817 28 C775 123 840 214 805 295 C786 339 792 390 824 432 L900 432 L900 28 Z"
            aria-hidden="true"
          />
          <path
            className={styles.frontierGlow}
            d="M817 28 C775 123 840 214 805 295 C786 339 792 390 824 432"
            aria-hidden="true"
          />
          <path
            className={styles.frontier}
            d="M817 28 C775 123 840 214 805 295 C786 339 792 390 824 432"
            aria-hidden="true"
          />
          <g className={styles.frontierLabel} aria-hidden="true">
            <text x="846" y="61">
              ORDER-23
            </text>
            <text x="846" y="77">
              FRONTIER
            </text>
            <circle cx="826" cy="67" r="3" />
          </g>

          <g className={styles.routes} aria-hidden="true">
            {historicalRoutePaths.map((path) => (
              <path d={path} key={path} />
            ))}
          </g>

          <g className={styles.routeProgress} aria-hidden="true">
            {routePaths.map((path, index) => {
              const arm = liveArms[index];
              const progress = arm
                ? Math.min(100, Math.max(0, arm.progress_percent))
                : 0;
              return (
                <path
                  d={path}
                  key={path}
                  pathLength="100"
                  style={{ strokeDasharray: `${progress} 100` }}
                />
              );
            })}
          </g>

          <g className={styles.basin} data-accent="cyan">
            <ellipse
              className={styles.radiusTwo}
              cx="212"
              cy="141"
              rx="109"
              ry="91"
            />
            <ellipse
              className={styles.radiusOne}
              cx="212"
              cy="141"
              rx="77"
              ry="63"
            />
            <BasinNodes nodes={alphaNodes} accent="cyan" />
            <g className={styles.basinLabel}>
              <text className={styles.basinName} x="116" y="61">
                GSDS COMPONENT
              </text>
              <text className={styles.basinCount} x="116" y="77">
                10 H · 7 HT
              </text>
            </g>
          </g>

          <g className={styles.basin} data-accent="violet">
            <ellipse
              className={styles.radiusTwo}
              cx="195"
              cy="327"
              rx="96"
              ry="79"
            />
            <ellipse
              className={styles.radiusOne}
              cx="195"
              cy="327"
              rx="64"
              ry="52"
            />
            <BasinNodes nodes={betaNodes} accent="violet" />
            <g className={styles.basinLabel}>
              <text className={styles.basinName} x="112" y="404">
                MENDELEY COMPONENT
              </text>
              <text className={styles.basinCount} x="112" y="420">
                4 H · 4 HT
              </text>
            </g>
          </g>

          <g className={styles.basin} data-accent="amber">
            <ellipse
              className={styles.harvestHalo}
              cx="445"
              cy="168"
              rx="72"
              ry="58"
            />
            <BasinNodes nodes={gammaNodes} accent="amber" />
            <g className={styles.basinLabel}>
              <text className={styles.basinName} x="382" y="232">
                HARVEST I
              </text>
              <text className={styles.basinCount} x="382" y="248">
                4 H · 4 HT
              </text>
            </g>
          </g>

          <g className={styles.basin} data-accent="rose">
            <ellipse
              className={styles.harvestHalo}
              cx="566"
              cy="318"
              rx="82"
              ry="64"
            />
            <BasinNodes nodes={deltaNodes} accent="rose" />
            <g className={styles.basinLabel}>
              <text className={styles.basinName} x="515" y="399">
                HARVEST II · NEW
              </text>
              <text className={styles.basinCount} x="515" y="415">
                6 H · 5 HT
              </text>
            </g>
          </g>

          <g className={styles.basin} data-accent="blue">
            <ellipse
              className={styles.harvestHalo}
              cx="648"
              cy="130"
              rx="76"
              ry="61"
            />
            <BasinNodes nodes={epsilonNodes} accent="blue" />
            <g className={styles.basinLabel}>
              <text className={styles.basinName} x="588" y="54">
                HARVEST III · NEW
              </text>
              <text className={styles.basinCount} x="588" y="70">
                6 H · 6 HT
              </text>
            </g>
          </g>

          <g className={styles.probes} aria-hidden="true">
            <g className={`${styles.probe} ${styles.probeOne}`}>
              <circle className={styles.probeWake} r="10" />
              <circle className={styles.probeCore} r="3.4" />
            </g>
            <g className={`${styles.probe} ${styles.probeTwo}`}>
              <circle className={styles.probeWake} r="10" />
              <circle className={styles.probeCore} r="3.4" />
            </g>
            <g className={`${styles.probe} ${styles.probeThree}`}>
              <circle className={styles.probeWake} r="10" />
              <circle className={styles.probeCore} r="3.4" />
            </g>
            <g className={`${styles.probe} ${styles.probeFour}`}>
              <circle className={styles.probeWake} r="10" />
              <circle className={styles.probeCore} r="3.4" />
            </g>
          </g>

          <g className={styles.probeLabel} aria-hidden="true">
            <path d="M658 111 L690 91" />
            <text x="696" y="88">
              {probeCount} NEW-BASIN TRAJECTORIES
            </text>
          </g>

          <g className={styles.axisLabels} aria-hidden="true">
            <text x="38" y="447">
              LOCAL STRUCTURE
            </text>
            <text x="714" y="447">
              STRICT-WIN DIRECTION →
            </text>
          </g>
        </svg>

        <div className={styles.legend} aria-label="Map legend">
          <span>
            <i className={styles.legendNode} aria-hidden="true" />
            H class
          </span>
          <span>
            <i className={styles.legendRadiusOne} aria-hidden="true" />
            exact r=1
          </span>
          <span>
            <i className={styles.legendRadiusTwo} aria-hidden="true" />
            exact r=2
          </span>
          <span>
            <i className={styles.legendHarvest} aria-hidden="true" />
            harvested basin
          </span>
          <span>
            <i className={styles.legendProbe} aria-hidden="true" />
            search probe
          </span>
        </div>
      </div>

      {liveArms.length > 0 ? (
        <div className={styles.liveRuns} aria-label="Live search telemetry">
          {liveArms.map((arm) => (
            <div className={styles.liveRun} key={arm.id} title={arm.engine}>
              <div className={styles.liveRunHead}>
                <strong>{arm.label}</strong>
                <span>{arm.best_ratio_percent.toFixed(2)}%</span>
              </div>
              <div
                className={styles.liveRunTrack}
                aria-label={`${arm.label}: ${arm.progress_percent.toFixed(1)}% of run time`}
              >
                <span
                  style={{
                    width: `${Math.min(100, Math.max(0, arm.progress_percent))}%`,
                  }}
                />
              </div>
              <div className={styles.liveRunMeta}>
                <span>{telemetryIsFresh ? arm.status : "snapshot"}</span>
                <span>
                  {formatDuration(arm.elapsed_seconds)} /{" "}
                  {formatDuration(arm.budget_seconds)}
                </span>
              </div>
            </div>
          ))}
        </div>
      ) : null}

      <dl className={styles.coverage} aria-label="Current search coverage">
        <div>
          <dt>Exact border r=1</dt>
          <dd>{radiusOneAssignments}</dd>
        </div>
        <div>
          <dt>Exact border r=2</dt>
          <dd>{radiusTwoAssignments}</dd>
        </div>
        <div>
          <dt>Basin trials</dt>
          <dd>{basinTrials}</dd>
        </div>
        <div>
          <dt>HT classes</dt>
          <dd>{knownHtClasses}</dd>
        </div>
      </dl>

      <footer className={styles.footer}>
        <p>
          <span className={styles.liveDot} aria-hidden="true" />
          {resolvedStatus}
        </p>
        <p>
          Schematic topology · route fill is run time · geometry is not to
          scale
        </p>
      </footer>

      <p className={styles.visuallyHidden}>
        Exact coverage has closed five seeded order-22 components containing
        30 H classes and 26 HT classes. Radius-one and radius-two border rings
        apply to the two published seed basins; three further basins were found
        by randomized harvesting and closed only under neutral one-entry moves.
        The map does not claim that no other components exist.
        Animated points represent active order-23 search probes approaching,
        but not crossing, the current determinant frontier. Filled route length
        represents elapsed run time, not mathematical distance.
      </p>
    </section>
  );
}
