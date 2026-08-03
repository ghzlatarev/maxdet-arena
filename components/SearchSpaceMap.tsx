"use client";

import { useEffect, useState, type CSSProperties } from "react";

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

type LegacySearchProgress = {
  schema_version?: 1;
  updated_at: string;
  frontier_absolute_determinant: string;
  arms: LiveSearchArm[];
};

type GramCampaignStatus =
  | "pending"
  | "active"
  | "complete"
  | "stopped"
  | "stale";

type GramCampaignHistoryPoint = {
  elapsed_seconds: number;
  best_absolute_determinant: string;
  best_ratio_percent: number;
  archive_size: number;
  sketch_discoveries: number;
  epochs_completed: number;
};

type GramCampaignArm = {
  id: string;
  label: string;
  engine: string;
  status: GramCampaignStatus;
  source_updated_at: string | null;
  elapsed_seconds: number;
  budget_seconds: number;
  progress_percent: number;
  best_absolute_determinant: string;
  best_ratio_percent: number;
  best_core_quotient: number;
  seed_basin_ids: string[];
  random_seed: number;
  kick_flips: number;
  epoch_moves: number;
  archive_size: number;
  archive_capacity: number;
  sketch_discoveries: number;
  epochs_completed: number;
  strict_target_states: number;
  history: GramCampaignHistoryPoint[];
};

type GramCampaignThreshold = {
  id: "archive-gate" | "near-frontier" | "frontier" | "strict";
  label: string;
  core_quotient: string;
  absolute_determinant: string;
};

type GramSeedBasin = {
  id: string;
  label: string;
  absolute_determinant: string;
  gram_basin_key_sha256: string;
};

type CoronalProjectionPoint = {
  arm_id: string;
  rank: number;
  orientation: number;
  core_quotient: string;
  absolute_determinant: string;
  det_m: string;
  kappa_numerator: string;
  kappa_denominator: string;
  kappa_decimal: number;
  pareto_front: boolean;
  identity?: "frontier" | "seed";
  seed_id?: string;
  seed_label?: string;
};

type CoronalProjection = {
  kind: "retained_exact_coronal";
  x_axis: "det_m";
  x_preference: "larger";
  y_axis: "kappa";
  y_preference: "smaller";
  points: CoronalProjectionPoint[];
};

type GramCampaignProgress = {
  schema_version: 2;
  updated_at: string;
  source_updated_at: string | null;
  campaign: {
    id: string;
    label: string;
    status: GramCampaignStatus;
  };
  frontier_absolute_determinant: string;
  strict_target_absolute_determinant: string;
  best_absolute_determinant: string;
  best_ratio_percent: number;
  thresholds: GramCampaignThreshold[];
  seed_basins: GramSeedBasin[];
  totals: {
    archive_cells: number;
    sketch_discoveries: number;
    epochs_completed: number;
    strict_target_states: number;
  };
  arms: GramCampaignArm[];
  coronal_projection?: CoronalProjection;
  claim_boundary: string;
};

const LIVE_TELEMETRY_MAX_AGE_MS = 15_000;
const GRAM_TELEMETRY_MAX_AGE_MS = 30_000;

type LiveSearchProgress = LegacySearchProgress | GramCampaignProgress;

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

function isRecord(value: unknown): value is Record<string, unknown> {
  return Boolean(value) && typeof value === "object";
}

function isFiniteNumber(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value);
}

function isExactInteger(value: unknown): value is string {
  return typeof value === "string" && /^\d+$/.test(value);
}

function isExactSignedInteger(value: unknown): value is string {
  return typeof value === "string" && /^-?\d+$/.test(value);
}

function isGramStatus(value: unknown): value is GramCampaignStatus {
  return (
    value === "pending" ||
    value === "active" ||
    value === "complete" ||
    value === "stopped" ||
    value === "stale"
  );
}

function isLegacySearchProgress(
  value: unknown,
): value is LegacySearchProgress {
  if (!isRecord(value)) {
    return false;
  }

  const candidate = value as Partial<LegacySearchProgress>;
  return (
    (candidate.schema_version === undefined ||
      candidate.schema_version === 1) &&
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

function isGramHistoryPoint(value: unknown): value is GramCampaignHistoryPoint {
  if (!isRecord(value)) return false;
  return (
    isFiniteNumber(value.elapsed_seconds) &&
    isExactInteger(value.best_absolute_determinant) &&
    isFiniteNumber(value.best_ratio_percent) &&
    isFiniteNumber(value.archive_size) &&
    isFiniteNumber(value.sketch_discoveries) &&
    isFiniteNumber(value.epochs_completed)
  );
}

function isGramCampaignArm(value: unknown): value is GramCampaignArm {
  if (!isRecord(value)) return false;
  return (
    typeof value.id === "string" &&
    typeof value.label === "string" &&
    typeof value.engine === "string" &&
    isGramStatus(value.status) &&
    (value.source_updated_at === null ||
      typeof value.source_updated_at === "string") &&
    isFiniteNumber(value.elapsed_seconds) &&
    isFiniteNumber(value.budget_seconds) &&
    isFiniteNumber(value.progress_percent) &&
    isExactInteger(value.best_absolute_determinant) &&
    isFiniteNumber(value.best_ratio_percent) &&
    isFiniteNumber(value.best_core_quotient) &&
    Array.isArray(value.seed_basin_ids) &&
    value.seed_basin_ids.every((item) => typeof item === "string") &&
    isFiniteNumber(value.random_seed) &&
    isFiniteNumber(value.kick_flips) &&
    isFiniteNumber(value.epoch_moves) &&
    isFiniteNumber(value.archive_size) &&
    isFiniteNumber(value.archive_capacity) &&
    isFiniteNumber(value.sketch_discoveries) &&
    isFiniteNumber(value.epochs_completed) &&
    isFiniteNumber(value.strict_target_states) &&
    Array.isArray(value.history) &&
    value.history.every(isGramHistoryPoint)
  );
}

function isCoronalProjectionPoint(
  value: unknown,
): value is CoronalProjectionPoint {
  if (!isRecord(value)) return false;
  return (
    typeof value.arm_id === "string" &&
    Number.isInteger(value.rank) &&
    Number.isInteger(value.orientation) &&
    isExactInteger(value.core_quotient) &&
    isExactInteger(value.absolute_determinant) &&
    isExactInteger(value.det_m) &&
    isExactSignedInteger(value.kappa_numerator) &&
    isExactInteger(value.kappa_denominator) &&
    BigInt(value.kappa_denominator) > 0n &&
    isFiniteNumber(value.kappa_decimal) &&
    typeof value.pareto_front === "boolean" &&
    (value.identity === undefined ||
      value.identity === "frontier" ||
      value.identity === "seed") &&
    (value.seed_id === undefined || typeof value.seed_id === "string") &&
    (value.seed_label === undefined || typeof value.seed_label === "string")
  );
}

function isCoronalProjection(value: unknown): value is CoronalProjection {
  if (!isRecord(value)) return false;
  return (
    value.kind === "retained_exact_coronal" &&
    value.x_axis === "det_m" &&
    value.x_preference === "larger" &&
    value.y_axis === "kappa" &&
    value.y_preference === "smaller" &&
    Array.isArray(value.points) &&
    value.points.length > 0 &&
    value.points.every(isCoronalProjectionPoint)
  );
}

function isGramCampaignProgress(
  value: unknown,
): value is GramCampaignProgress {
  if (!isRecord(value) || value.schema_version !== 2) return false;
  if (!isRecord(value.campaign) || !isRecord(value.totals)) return false;
  return (
    typeof value.updated_at === "string" &&
    (value.source_updated_at === null ||
      typeof value.source_updated_at === "string") &&
    typeof value.campaign.id === "string" &&
    typeof value.campaign.label === "string" &&
    isGramStatus(value.campaign.status) &&
    isExactInteger(value.frontier_absolute_determinant) &&
    isExactInteger(value.strict_target_absolute_determinant) &&
    isExactInteger(value.best_absolute_determinant) &&
    isFiniteNumber(value.best_ratio_percent) &&
    Array.isArray(value.thresholds) &&
    value.thresholds.every(
      (threshold) =>
        isRecord(threshold) &&
        (threshold.id === "archive-gate" ||
          threshold.id === "near-frontier" ||
          threshold.id === "frontier" ||
          threshold.id === "strict") &&
        typeof threshold.label === "string" &&
        isExactInteger(threshold.core_quotient) &&
        isExactInteger(threshold.absolute_determinant),
    ) &&
    Array.isArray(value.seed_basins) &&
    value.seed_basins.every(
      (seed) =>
        isRecord(seed) &&
        typeof seed.id === "string" &&
        typeof seed.label === "string" &&
        isExactInteger(seed.absolute_determinant) &&
        typeof seed.gram_basin_key_sha256 === "string",
    ) &&
    isFiniteNumber(value.totals.archive_cells) &&
    isFiniteNumber(value.totals.sketch_discoveries) &&
    isFiniteNumber(value.totals.epochs_completed) &&
    isFiniteNumber(value.totals.strict_target_states) &&
    Array.isArray(value.arms) &&
    value.arms.every(isGramCampaignArm) &&
    (value.coronal_projection === undefined ||
      isCoronalProjection(value.coronal_projection)) &&
    typeof value.claim_boundary === "string"
  );
}

function isLiveSearchProgress(value: unknown): value is LiveSearchProgress {
  return isGramCampaignProgress(value) || isLegacySearchProgress(value);
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

const GRAM_ARM_COLORS = [
  "#55e6d2",
  "#9e87ff",
  "#ffbd62",
  "#ff7895",
  "#66a7ff",
  "#8ee86f",
] as const;

const gramThresholdY: Record<GramCampaignThreshold["id"], number> = {
  "archive-gate": 396,
  "near-frontier": 286,
  frontier: 137,
  strict: 67,
};

const GRAM_GRAPH_LEFT = 164;
const GRAM_GRAPH_RIGHT = 884;
const GRAM_GRAPH_WIDTH = GRAM_GRAPH_RIGHT - GRAM_GRAPH_LEFT;

function formatExactInteger(value: string | number): string {
  return BigInt(value).toLocaleString("en-US");
}

function clamp(value: number, minimum: number, maximum: number): number {
  return Math.min(maximum, Math.max(minimum, value));
}

function gramScoreY(
  score: string,
  thresholds: GramCampaignThreshold[],
): number {
  const scale = (
    ["archive-gate", "near-frontier", "frontier", "strict"] as const
  )
    .map((id) => {
      const threshold = thresholds.find((item) => item.id === id);
      return threshold
        ? {
            score: Number(threshold.absolute_determinant),
            y: gramThresholdY[id],
          }
        : null;
    })
    .filter(
      (
        item,
      ): item is {
        score: number;
        y: number;
      } => item !== null,
    );
  const value = Number(score);
  if (scale.length < 2 || !Number.isFinite(value)) return 396;
  if (value <= scale[0].score) return scale[0].y;
  for (let index = 1; index < scale.length; index += 1) {
    const lower = scale[index - 1];
    const upper = scale[index];
    if (value <= upper.score) {
      const span = upper.score - lower.score;
      const fraction = span > 0 ? (value - lower.score) / span : 1;
      return lower.y + fraction * (upper.y - lower.y);
    }
  }
  return scale[scale.length - 1].y;
}

function gramArmStatus(
  arm: GramCampaignArm,
  telemetryClock: number,
): GramCampaignStatus {
  if (arm.status !== "active") return arm.status;
  if (!arm.source_updated_at) return "stale";
  const age = telemetryClock - Date.parse(arm.source_updated_at);
  return Number.isFinite(age) &&
    age >= -5_000 &&
    age <= GRAM_TELEMETRY_MAX_AGE_MS
    ? "active"
    : "stale";
}

function aggregateGramStatus(statuses: GramCampaignStatus[]): string {
  if (statuses.length === 0) return "Waiting for search arms";
  if (statuses.every((status) => status === "complete")) {
    return `${statuses.length} search arms complete`;
  }
  const active = statuses.filter((status) => status === "active").length;
  if (active > 0) return `${active} search arms active`;
  if (statuses.includes("stale")) return "Search telemetry stale";
  if (statuses.includes("stopped")) return "Search stopped";
  return "Waiting for search arms";
}

const CORONAL_PLOT_LEFT = 94;
const CORONAL_PLOT_RIGHT = 884;
const CORONAL_PLOT_TOP = 48;
const CORONAL_PLOT_BOTTOM = 424;
const CORONAL_PLOT_WIDTH = CORONAL_PLOT_RIGHT - CORONAL_PLOT_LEFT;
const CORONAL_PLOT_HEIGHT = CORONAL_PLOT_BOTTOM - CORONAL_PLOT_TOP;

function formatDetMScale(value: number): string {
  const scaled = value / 1e15;
  return scaled >= 100 ? scaled.toFixed(0) : scaled.toFixed(1);
}

function coronalArmLabel(arm: GramCampaignArm): string {
  return arm.id
    .replace(/^arm-/, "")
    .replace(/-k(\d+)$/, " · k=$1")
    .replaceAll("-", " ");
}

function CoronalCampaignMap({
  className,
  progress,
  telemetryClock,
}: {
  className?: string;
  progress: GramCampaignProgress & {
    coronal_projection: CoronalProjection;
  };
  telemetryClock: number;
}) {
  const rootClassName = [
    styles.root,
    styles.gramRoot,
    styles.coronalRoot,
    className,
  ]
    .filter(Boolean)
    .join(" ");
  const statuses = progress.arms.map((arm) =>
    gramArmStatus(arm, telemetryClock),
  );
  const statusLabel = aggregateGramStatus(statuses);
  const points = progress.coronal_projection.points;
  const retainedMatrices = new Set(
    points.map((point) => `${point.arm_id}:${point.rank}`),
  ).size;
  const paretoMatrices = new Set(
    points
      .filter((point) => point.pareto_front)
      .map((point) => `${point.arm_id}:${point.rank}`),
  ).size;
  const armLookup = new Map(
    progress.arms.map((arm, index) => [
      arm.id,
      {
        label: coronalArmLabel(arm),
        color: GRAM_ARM_COLORS[index % GRAM_ARM_COLORS.length],
      },
    ]),
  );

  const detMValues = points.map((point) => Number(point.det_m));
  const kappaValues = points.map((point) => point.kappa_decimal);
  const rawMinDetM = Math.min(...detMValues);
  const rawMaxDetM = Math.max(...detMValues);
  const rawMinKappa = Math.min(...kappaValues);
  const rawMaxKappa = Math.max(...kappaValues);
  const detMSpan = Math.max(1, rawMaxDetM - rawMinDetM);
  const kappaSpan = Math.max(0.01, rawMaxKappa - rawMinKappa);
  const minDetM = rawMinDetM - detMSpan * 0.055;
  const maxDetM = rawMaxDetM + detMSpan * 0.055;
  const minKappa = rawMinKappa - kappaSpan * 0.11;
  const maxKappa = rawMaxKappa + kappaSpan * 0.11;
  const plotX = (value: number) =>
    CORONAL_PLOT_LEFT +
    ((value - minDetM) / (maxDetM - minDetM)) * CORONAL_PLOT_WIDTH;
  // A smaller kappa is desirable and therefore appears higher in the plane.
  const plotY = (value: number) =>
    CORONAL_PLOT_TOP +
    ((value - minKappa) / (maxKappa - minKappa)) * CORONAL_PLOT_HEIGHT;
  const xTicks = Array.from(
    { length: 5 },
    (_, index) => minDetM + (index * (maxDetM - minDetM)) / 4,
  );
  const yTicks = Array.from(
    { length: 5 },
    (_, index) => minKappa + (index * (maxKappa - minKappa)) / 4,
  );

  const locationTotals = new Map<string, number>();
  for (const point of points) {
    const key = `${point.det_m}:${point.kappa_numerator}/${point.kappa_denominator}`;
    locationTotals.set(key, (locationTotals.get(key) ?? 0) + 1);
  }
  const locationIndices = new Map<string, number>();
  const plottedPoints = points.map((point) => {
    const key = `${point.det_m}:${point.kappa_numerator}/${point.kappa_denominator}`;
    const duplicateCount = locationTotals.get(key) ?? 1;
    const duplicateIndex = locationIndices.get(key) ?? 0;
    locationIndices.set(key, duplicateIndex + 1);
    const angle = (2 * Math.PI * duplicateIndex) / duplicateCount;
    const offset = duplicateCount > 1 ? Math.min(6, 2.2 + duplicateCount) : 0;
    const x = plotX(Number(point.det_m));
    const y = plotY(point.kappa_decimal);
    return {
      ...point,
      x,
      y,
      displayX: x + Math.cos(angle) * offset,
      displayY: y + Math.sin(angle) * offset,
    };
  });

  const paretoByCoordinate = new Map<
    string,
    (typeof plottedPoints)[number]
  >();
  for (const point of plottedPoints) {
    if (!point.pareto_front) continue;
    const key = `${point.det_m}:${point.kappa_numerator}/${point.kappa_denominator}`;
    if (!paretoByCoordinate.has(key)) {
      paretoByCoordinate.set(key, point);
    }
  }
  const paretoPoints = [...paretoByCoordinate.values()].sort(
    (left, right) => Number(left.det_m) - Number(right.det_m),
  );
  const firstFrontierIndex = plottedPoints.findIndex(
    (point) => point.identity === "frontier",
  );

  return (
    <section
      className={rootClassName}
      aria-label="Exact coronal projection of retained MaxDet matrices"
    >
      <header className={styles.header}>
        <div>
          <p className={styles.eyebrow}>
            Coronal plane · retained exact matrices
          </p>
          <h3>The archive, viewed in a different geometry.</h3>
        </div>
        <dl
          className={`${styles.metrics} ${styles.gramMetrics}`}
          aria-label="Exact coronal campaign totals"
        >
          <div>
            <dt>Seed basins</dt>
            <dd>{progress.seed_basins.length}</dd>
          </div>
          <div>
            <dt>Retained matrices</dt>
            <dd>{retainedMatrices}</dd>
          </div>
          <div className={styles.gramBestMetric}>
            <dt>Exact best</dt>
            <dd title={formatExactInteger(progress.best_absolute_determinant)}>
              {formatExactInteger(progress.best_absolute_determinant)}
            </dd>
          </div>
        </dl>
      </header>

      <div className={`${styles.mapFrame} ${styles.coronalMapFrame}`}>
        <svg
          className={`${styles.map} ${styles.coronalMap}`}
          viewBox="0 0 920 480"
          role="img"
          aria-label={`${retainedMatrices} retained exact matrices projected by det M and kappa; larger det M is right, smaller kappa is up`}
          preserveAspectRatio="xMidYMid meet"
        >
          <title>Exact coronal-Pareto archive projection</title>
          <desc>
            This is a two-dimensional projection of retained exact matrices,
            not coverage of the full search space. Horizontal position is exact
            det M, increasing to the right. Vertical position is exact kappa,
            with smaller values higher. Colors identify search arms. A bright
            line joins the nondominated points in these two plotted
            coordinates, and the published frontier matrix is ringed.
          </desc>

          <rect
            className={styles.coronalPlotField}
            x={CORONAL_PLOT_LEFT}
            y={CORONAL_PLOT_TOP}
            width={CORONAL_PLOT_WIDTH}
            height={CORONAL_PLOT_HEIGHT}
            rx="14"
            aria-hidden="true"
          />

          <g className={styles.coronalGrid} aria-hidden="true">
            {xTicks.map((tick, index) => {
              const x = plotX(tick);
              return (
                <g key={`coronal-x-${index}`}>
                  <line
                    x1={x}
                    x2={x}
                    y1={CORONAL_PLOT_TOP}
                    y2={CORONAL_PLOT_BOTTOM}
                  />
                  <text
                    textAnchor="middle"
                    x={x}
                    y={CORONAL_PLOT_BOTTOM + 22}
                  >
                    {formatDetMScale(tick)}
                  </text>
                </g>
              );
            })}
            {yTicks.map((tick, index) => {
              const y = plotY(tick);
              return (
                <g key={`coronal-y-${index}`}>
                  <line
                    x1={CORONAL_PLOT_LEFT}
                    x2={CORONAL_PLOT_RIGHT}
                    y1={y}
                    y2={y}
                  />
                  <text
                    textAnchor="end"
                    x={CORONAL_PLOT_LEFT - 12}
                    y={y + 3}
                  >
                    {tick.toFixed(3)}
                  </text>
                </g>
              );
            })}
          </g>

          {paretoPoints.length > 1 ? (
            <polyline
              className={styles.coronalParetoLine}
              points={paretoPoints
                .map((point) => `${point.x},${point.y}`)
                .join(" ")}
              aria-hidden="true"
            />
          ) : null}

          <g role="list" aria-label="Retained exact matrix projections">
            {plottedPoints.map((point, index) => {
              const armDetails = armLookup.get(point.arm_id);
              const armLabel = armDetails?.label ?? point.arm_id;
              const armColor =
                armDetails?.color ?? GRAM_ARM_COLORS[index % GRAM_ARM_COLORS.length];
              const pointStyle = {
                "--arm-color": armColor,
              } as CSSProperties;
              const identityLabel =
                point.identity === "frontier"
                  ? "published frontier"
                  : point.identity === "seed"
                    ? point.seed_label ?? "campaign seed"
                    : null;
              const accessibleLabel = `${armLabel}, archive rank ${
                point.rank
              }, absolute determinant ${formatExactInteger(
                point.absolute_determinant,
              )}, det M ${formatExactInteger(point.det_m)}, kappa ${
                point.kappa_numerator
              } over ${point.kappa_denominator}${
                point.pareto_front ? ", on the plotted Pareto front" : ""
              }${identityLabel ? `, ${identityLabel}` : ""}`;
              return (
                <g
                  className={styles.coronalPointGroup}
                  data-frontier={
                    point.identity === "frontier" ? "true" : undefined
                  }
                  data-pareto={point.pareto_front ? "true" : undefined}
                  data-seed={point.identity === "seed" ? "true" : undefined}
                  key={`${point.arm_id}-${point.rank}-${point.orientation}`}
                  role="listitem"
                  aria-label={accessibleLabel}
                  style={pointStyle}
                  tabIndex={0}
                >
                  <title>{accessibleLabel}</title>
                  {point.pareto_front ? (
                    <circle
                      className={styles.coronalParetoHalo}
                      cx={point.displayX}
                      cy={point.displayY}
                      r="7.5"
                    />
                  ) : null}
                  {point.identity === "seed" ? (
                    <circle
                      className={styles.coronalSeedRing}
                      cx={point.displayX}
                      cy={point.displayY}
                      r="6.2"
                    />
                  ) : null}
                  {point.identity === "frontier" ? (
                    <circle
                      className={styles.coronalFrontierRing}
                      cx={point.displayX}
                      cy={point.displayY}
                      r="9.5"
                    />
                  ) : null}
                  <circle
                    className={styles.coronalPoint}
                    cx={point.displayX}
                    cy={point.displayY}
                    r={point.pareto_front ? 4.2 : 3.2}
                  />
                  {index === firstFrontierIndex ? (
                    <g className={styles.coronalFrontierLabel}>
                      <line
                        x1={point.displayX + 7}
                        x2={point.displayX + 35}
                        y1={point.displayY - 7}
                        y2={point.displayY - 25}
                      />
                      <text
                        x={point.displayX + 39}
                        y={point.displayY - 27}
                      >
                        FRONTIER MATRIX
                      </text>
                    </g>
                  ) : null}
                </g>
              );
            })}
          </g>

          <g className={styles.coronalAxes} aria-hidden="true">
            <text
              textAnchor="middle"
              x={(CORONAL_PLOT_LEFT + CORONAL_PLOT_RIGHT) / 2}
              y="474"
            >
              det(M) × 10¹⁵ · LARGER →
            </text>
            <text
              textAnchor="middle"
              transform={`translate(20 ${
                (CORONAL_PLOT_TOP + CORONAL_PLOT_BOTTOM) / 2
              }) rotate(-90)`}
            >
              κ · SMALLER ↑
            </text>
            <text x={CORONAL_PLOT_LEFT + 12} y={CORONAL_PLOT_TOP + 20}>
              DESIRABLE ↑
            </text>
            <text
              textAnchor="end"
              x={CORONAL_PLOT_RIGHT - 12}
              y={CORONAL_PLOT_BOTTOM - 12}
            >
              RETAINED ARCHIVE PROJECTION
            </text>
          </g>
        </svg>

        <div
          className={styles.coronalLegend}
          aria-label="Coronal projection legend"
        >
          {progress.arms.map((arm, index) => {
            const legendStyle = {
              "--arm-color":
                GRAM_ARM_COLORS[index % GRAM_ARM_COLORS.length],
            } as CSSProperties;
            return (
              <span key={arm.id} style={legendStyle}>
                <i className={styles.coronalLegendArm} aria-hidden="true" />
                {coronalArmLabel(arm)}
              </span>
            );
          })}
          <span>
            <i className={styles.coronalLegendSeed} aria-hidden="true" />
            campaign seed
          </span>
          <span>
            <i className={styles.coronalLegendPareto} aria-hidden="true" />
            plotted Pareto front
          </span>
          <span>
            <i className={styles.coronalLegendFrontier} aria-hidden="true" />
            frontier matrix
          </span>
        </div>
      </div>

      <div
        className={styles.coronalArms}
        aria-label="Coronal-Pareto campaign arms"
      >
        {progress.arms.map((arm, index) => {
          const armStyle = {
            "--arm-color":
              GRAM_ARM_COLORS[index % GRAM_ARM_COLORS.length],
          } as CSSProperties;
          return (
            <article className={styles.coronalArm} key={arm.id} style={armStyle}>
              <div>
                <strong title={arm.id}>{coronalArmLabel(arm)}</strong>
                <span data-status={statuses[index]}>{statuses[index]}</span>
              </div>
              <p>
                {arm.archive_size} retained · best{" "}
                {formatExactInteger(arm.best_absolute_determinant)}
              </p>
            </article>
          );
        })}
      </div>

      <dl className={styles.coverage} aria-label="Coronal projection totals">
        <div>
          <dt>Retained matrices</dt>
          <dd>{retainedMatrices}</dd>
        </div>
        <div>
          <dt>Exact orientations</dt>
          <dd>{points.length}</dd>
        </div>
        <div>
          <dt>2D Pareto matrices</dt>
          <dd>{paretoMatrices}</dd>
        </div>
        <div>
          <dt>Strict target states</dt>
          <dd>{formatExactInteger(progress.totals.strict_target_states)}</dd>
        </div>
      </dl>

      <footer className={styles.footer}>
        <p>
          <span className={styles.liveDot} aria-hidden="true" />
          {statusLabel}
        </p>
        <p>Projection of retained exact matrices · not global coverage</p>
      </footer>

      <p className={styles.gramClaim}>{progress.claim_boundary}</p>
    </section>
  );
}

function GramCampaignMap({
  className,
  progress,
  telemetryClock,
}: {
  className?: string;
  progress: GramCampaignProgress;
  telemetryClock: number;
}) {
  const rootClassName = [styles.root, styles.gramRoot, className]
    .filter(Boolean)
    .join(" ");
  const statuses = progress.arms.map((arm) =>
    gramArmStatus(arm, telemetryClock),
  );
  const statusLabel = aggregateGramStatus(statuses);
  const thresholdLookup = new Map(
    progress.thresholds.map((threshold) => [threshold.id, threshold]),
  );

  if (
    progress.coronal_projection &&
    progress.coronal_projection.points.length > 0
  ) {
    return (
      <CoronalCampaignMap
        className={className}
        progress={
          progress as GramCampaignProgress & {
            coronal_projection: CoronalProjection;
          }
        }
        telemetryClock={telemetryClock}
      />
    );
  }

  return (
    <section
      className={rootClassName}
      aria-label="Gram-sketch basin-hopper campaign telemetry"
    >
      <header className={styles.header}>
        <div>
          <p className={styles.eyebrow}>
            Gram-basin campaign · exact telemetry
          </p>
          <h3>Seed basins fan into a live sketch archive.</h3>
        </div>
        <dl
          className={`${styles.metrics} ${styles.gramMetrics}`}
          aria-label="Gram campaign totals"
        >
          <div>
            <dt>Seed basins</dt>
            <dd>{progress.seed_basins.length}</dd>
          </div>
          <div>
            <dt>Archive cells</dt>
            <dd>{formatExactInteger(progress.totals.archive_cells)}</dd>
          </div>
          <div className={styles.gramBestMetric}>
            <dt>Exact best</dt>
            <dd title={formatExactInteger(progress.best_absolute_determinant)}>
              {formatExactInteger(progress.best_absolute_determinant)}
            </dd>
          </div>
        </dl>
      </header>

      <div className={`${styles.mapFrame} ${styles.gramMapFrame}`}>
        <svg
          className={`${styles.map} ${styles.gramMap}`}
          viewBox="0 0 920 460"
          role="img"
          aria-label={`${progress.seed_basins.length} exact local seed basins and ${progress.arms.length} Gram-sketch search arms plotted over elapsed time and schematic score bands`}
          preserveAspectRatio="xMidYMid meet"
        >
          <title>Gram-sketch basin-hopper campaign</title>
          <desc>
            Exact local seed basins begin at the left. Each colored arm moves
            across elapsed search time. Circles mark increases in the
            noncanonical Gram-sketch discovery count. Vertical position uses
            schematic score bands so the frontier and first strict score remain
            visually distinct.
          </desc>

          <g className={styles.grid} aria-hidden="true">
            {Array.from({ length: 13 }, (_, index) => (
              <line
                key={`gram-vertical-${index}`}
                x1={GRAM_GRAPH_LEFT + index * (GRAM_GRAPH_WIDTH / 12)}
                x2={GRAM_GRAPH_LEFT + index * (GRAM_GRAPH_WIDTH / 12)}
                y1="42"
                y2="421"
              />
            ))}
          </g>

          <g className={styles.gramBands} aria-hidden="true">
            <rect
              className={styles.gramStrictBand}
              x={GRAM_GRAPH_LEFT}
              y="43"
              width={GRAM_GRAPH_WIDTH}
              height={gramThresholdY.frontier - 43}
            />
            <rect
              className={styles.gramNearBand}
              x={GRAM_GRAPH_LEFT}
              y={gramThresholdY.frontier}
              width={GRAM_GRAPH_WIDTH}
              height={gramThresholdY["near-frontier"] - gramThresholdY.frontier}
            />
            {progress.thresholds.map((threshold) => (
              <g key={threshold.id} data-threshold={threshold.id}>
                <line
                  x1={GRAM_GRAPH_LEFT}
                  x2={GRAM_GRAPH_RIGHT}
                  y1={gramThresholdY[threshold.id]}
                  y2={gramThresholdY[threshold.id]}
                />
                <text
                  className={styles.gramThresholdLabel}
                  x={GRAM_GRAPH_LEFT + 10}
                  y={gramThresholdY[threshold.id] - 8}
                >
                  {threshold.label}
                </text>
                <text
                  className={styles.gramThresholdValue}
                  textAnchor="end"
                  x={GRAM_GRAPH_RIGHT - 8}
                  y={gramThresholdY[threshold.id] - 8}
                >
                  {formatExactInteger(threshold.absolute_determinant)}
                </text>
              </g>
            ))}
          </g>

          <g className={styles.gramSeedCloud}>
            <text className={styles.gramSeedHeading} x="33" y="33">
              EXACT LOCAL
            </text>
            <text className={styles.gramSeedHeading} x="33" y="47">
              SEED BASINS
            </text>
            {progress.seed_basins.map((seed, index) => {
              const x = 62 + (index % 3) * 25;
              const y =
                gramScoreY(seed.absolute_determinant, progress.thresholds) +
                (index % 2 === 0 ? -4 : 4);
              return (
                <g key={seed.id}>
                  <circle className={styles.gramSeedHalo} cx={x} cy={y} r="9" />
                  <circle className={styles.gramSeedCore} cx={x} cy={y} r="3.5">
                    <title>
                      {seed.label} ·{" "}
                      {formatExactInteger(seed.absolute_determinant)} · basin{" "}
                      {seed.gram_basin_key_sha256.slice(0, 12)}
                    </title>
                  </circle>
                </g>
              );
            })}
          </g>

          <line
            className={styles.gramSweep}
            x1={GRAM_GRAPH_LEFT}
            x2={GRAM_GRAPH_LEFT}
            y1="43"
            y2="421"
            aria-hidden="true"
          />

          <g className={styles.gramTrails}>
            {progress.arms.map((arm, armIndex) => {
              const history =
                arm.history.length > 0
                  ? arm.history
                  : [
                      {
                        elapsed_seconds: arm.elapsed_seconds,
                        best_absolute_determinant:
                          arm.best_absolute_determinant,
                        best_ratio_percent: arm.best_ratio_percent,
                        archive_size: arm.archive_size,
                        sketch_discoveries: arm.sketch_discoveries,
                        epochs_completed: arm.epochs_completed,
                      },
                    ];
              const laneOffset =
                (armIndex - (progress.arms.length - 1) / 2) * 6;
              const points = history.map((point) => {
                const fraction =
                  arm.budget_seconds > 0
                    ? clamp(
                        point.elapsed_seconds / arm.budget_seconds,
                        0,
                        1,
                      )
                    : 0;
                return {
                  ...point,
                  x: GRAM_GRAPH_LEFT + fraction * GRAM_GRAPH_WIDTH,
                  y:
                    gramScoreY(
                      point.best_absolute_determinant,
                      progress.thresholds,
                    ) + laneOffset,
                };
              });
              const firstPoint = points[0];
              const finalPoint = points[points.length - 1];
              const armColor = GRAM_ARM_COLORS[
                armIndex % GRAM_ARM_COLORS.length
              ];
              const armStyle = {
                "--arm-color": armColor,
              } as CSSProperties;
              const resolvedStatus = statuses[armIndex];
              return (
                <g key={arm.id} style={armStyle}>
                  <path
                    className={styles.gramSeedLink}
                    d={`M126 ${firstPoint.y} L${firstPoint.x} ${firstPoint.y}`}
                  />
                  <polyline
                    className={styles.gramTrail}
                    points={points
                      .map((point) => `${point.x},${point.y}`)
                      .join(" ")}
                  />
                  {points.map((point, pointIndex) => {
                    const previous =
                      pointIndex > 0
                        ? points[pointIndex - 1].sketch_discoveries
                        : 0;
                    const discoveries = Math.max(
                      0,
                      point.sketch_discoveries - previous,
                    );
                    if (discoveries === 0) return null;
                    const radius = Math.min(
                      5.4,
                      1.6 + Math.log2(discoveries + 1) * 0.48,
                    );
                    return (
                      <circle
                        className={styles.gramDiscovery}
                        cx={point.x}
                        cy={point.y}
                        data-recent={
                          pointIndex >= points.length - 6 ? "true" : undefined
                        }
                        key={`${arm.id}-${point.elapsed_seconds}`}
                        r={radius}
                      >
                        <title>
                          {arm.label} · {formatDuration(point.elapsed_seconds)} ·{" "}
                          {formatExactInteger(point.sketch_discoveries)}{" "}
                          cumulative Gram-sketch discoveries
                        </title>
                      </circle>
                    );
                  })}
                  <g
                    className={styles.gramHead}
                    data-active={
                      resolvedStatus === "active" ? "true" : undefined
                    }
                    transform={`translate(${finalPoint.x} ${finalPoint.y})`}
                  >
                    <circle className={styles.gramHeadWake} r="10" />
                    <circle className={styles.gramHeadCore} r="3.7" />
                  </g>
                  <text
                    className={styles.gramArmLabel}
                    textAnchor="end"
                    x={GRAM_GRAPH_RIGHT - 7}
                    y={finalPoint.y + 18}
                  >
                    {arm.label}
                  </text>
                </g>
              );
            })}
          </g>

          <g className={styles.axisLabels} aria-hidden="true">
            <text x={GRAM_GRAPH_LEFT} y="447">
              ELAPSED SEARCH TIME →
            </text>
            <text textAnchor="end" x={GRAM_GRAPH_RIGHT} y="447">
              SCORE BANDS ARE SCHEMATIC
            </text>
          </g>
        </svg>

        <div className={styles.gramLegend} aria-label="Campaign map legend">
          <span>
            <i className={styles.gramLegendSeed} aria-hidden="true" />
            exact seed basin
          </span>
          <span>
            <i className={styles.gramLegendDiscovery} aria-hidden="true" />
            sketch discovery
          </span>
          <span>
            <i className={styles.gramLegendTrail} aria-hidden="true" />
            search arm
          </span>
          <span>
            <i className={styles.gramLegendThreshold} aria-hidden="true" />
            score threshold
          </span>
        </div>
      </div>

      <div
        className={styles.gramLiveRuns}
        aria-label="Gram-sketch search arms"
      >
        {progress.arms.map((arm, index) => {
          const armStyle = {
            "--arm-color":
              GRAM_ARM_COLORS[index % GRAM_ARM_COLORS.length],
          } as CSSProperties;
          return (
            <article
              className={styles.gramLiveRun}
              key={arm.id}
              style={armStyle}
              title={arm.engine}
            >
              <header className={styles.gramLiveRunHead}>
                <strong>{arm.label}</strong>
                <span data-status={statuses[index]}>{statuses[index]}</span>
              </header>
              <p className={styles.gramRunScore}>
                {formatExactInteger(arm.best_absolute_determinant)}
              </p>
              <div
                className={styles.gramLiveRunTrack}
                aria-label={`${arm.label}: ${arm.progress_percent.toFixed(1)}% of run time`}
              >
                <span
                  style={{
                    width: `${clamp(arm.progress_percent, 0, 100)}%`,
                  }}
                />
              </div>
              <dl className={styles.gramRunMetrics}>
                <div>
                  <dt>Archive</dt>
                  <dd>
                    {arm.archive_size}/{arm.archive_capacity}
                  </dd>
                </div>
                <div>
                  <dt>Sketches</dt>
                  <dd>{formatExactInteger(arm.sketch_discoveries)}</dd>
                </div>
                <div>
                  <dt>Epochs</dt>
                  <dd>{formatExactInteger(arm.epochs_completed)}</dd>
                </div>
                <div>
                  <dt>Time</dt>
                  <dd>
                    {formatDuration(arm.elapsed_seconds)} /{" "}
                    {formatDuration(arm.budget_seconds)}
                  </dd>
                </div>
              </dl>
            </article>
          );
        })}
      </div>

      <dl className={styles.coverage} aria-label="Gram campaign totals">
        <div>
          <dt>Exact seed basins</dt>
          <dd>{progress.seed_basins.length}</dd>
        </div>
        <div>
          <dt>Arm-local sketches</dt>
          <dd>{formatExactInteger(progress.totals.sketch_discoveries)}</dd>
        </div>
        <div>
          <dt>Epochs completed</dt>
          <dd>{formatExactInteger(progress.totals.epochs_completed)}</dd>
        </div>
        <div>
          <dt>Strict target states</dt>
          <dd>{formatExactInteger(progress.totals.strict_target_states)}</dd>
        </div>
      </dl>

      <footer className={styles.footer}>
        <p>
          <span className={styles.liveDot} aria-hidden="true" />
          {statusLabel}
        </p>
        <p>Noncanonical Gram sketches · arm counts may overlap</p>
      </footer>

      <p className={styles.gramClaim}>{progress.claim_boundary}</p>
      <p className={styles.visuallyHidden}>
        The map distinguishes exactly classified local seed basins from
        noncanonical Gram-sketch search descriptors. Archive and discovery
        totals are summed across arms and may overlap. A different sketch
        proves a different signed Gram orbit, while a matching sketch does not
        prove equivalence.
        {thresholdLookup.get("frontier")
          ? ` The exact frontier is ${formatExactInteger(
              thresholdLookup.get("frontier")!.absolute_determinant,
            )}.`
          : ""}
      </p>
    </section>
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
    const arenaRoute = "/problems/maxdet";
    const arenaRouteIndex = window.location.pathname.indexOf(arenaRoute);
    const basePath =
      arenaRouteIndex >= 0
        ? window.location.pathname.slice(0, arenaRouteIndex)
        : "";
    const progressUrl = new URL(
      `${basePath}/search-progress.json`,
      window.location.origin,
    );

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

  if (liveProgress?.schema_version === 2) {
    return (
      <GramCampaignMap
        className={className}
        progress={liveProgress}
        telemetryClock={telemetryClock}
      />
    );
  }

  const legacyProgress = liveProgress;
  const liveArms =
    legacyProgress?.arms.slice(0, routePaths.length) ?? [];
  const telemetryAge = legacyProgress
    ? telemetryClock - Date.parse(legacyProgress.updated_at)
    : Number.POSITIVE_INFINITY;
  const telemetryIsFresh =
    Number.isFinite(telemetryAge) &&
    telemetryAge >= 0 &&
    telemetryAge <= LIVE_TELEMETRY_MAX_AGE_MS;
  const liveBestRatio =
    liveArms.length === 0
      ? bestRatio
      : `${Math.max(
          ...liveArms.map((arm) => arm.best_ratio_percent),
        ).toFixed(2)}%`;
  const activeSearches = telemetryIsFresh
    ? liveArms.filter((arm) => arm.status === "active").length
    : 0;
  const allSearchesComplete =
    liveArms.length > 0 &&
    liveArms.every((arm) => arm.status === "complete");
  const resolvedStatus =
    liveArms.length > 0
      ? allSearchesComplete
        ? `${liveArms.length} searches complete`
        : !telemetryIsFresh
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
                <span>
                  {arm.status === "complete" || telemetryIsFresh
                    ? arm.status
                    : "snapshot"}
                </span>
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
