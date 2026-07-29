type Audit = {
  label: string;
  assignments: string;
  detail: string;
};

type CampaignExtension = {
  id: string;
  kicker: string;
  headline: {
    value: string;
    label: string;
  };
  metrics: {
    label: string;
    value: string;
  }[];
  summary: string;
  evidence: string;
};

type PortalCampaign = {
  method: string;
  local_h_ht_classes_before: number;
  new_local_h_ht_classes: number;
  local_h_ht_classes_after: number;
  substantive_frontier_ties: number;
  strict_wins: number;
  exact_radius_3_audits: {
    id: string;
    label: string;
    assignments: string;
    frontier_ties: number;
    strict_wins: number;
    status: string;
  }[];
  exact_connector_cubes: {
    count: number;
    assignments: string;
    frontier_ties: number;
    strict_wins: number;
    result: string;
  };
  claim_boundary: string;
};

type Snapshot = {
  updated: string;
  status: string;
  portal_campaign: PortalCampaign;
  exact_audits: Audit[];
  radius_4_screen: Audit;
  h24_seed_sweep: {
    pivots: string;
    classes: number;
    frontier_ties: number;
    strict_wins: number;
  };
  neutral_frontier: {
    verified_raw_matrices: number;
    new_raw_matrices: number;
    exact_cycles: number;
    states_per_cycle: number;
    neutral_edges: number;
    max_hamming_from_baseline: number;
    strict_wins: number;
    claim_boundary: string;
  };
  qubo_trust_region: {
    models: string;
    exact_pair_coefficients: string;
    exact_block_evaluations: string;
    frontier_ties: number;
    new_local_h_classes: number;
    local_h_classes_after: number;
    strict_wins: number;
    claim_boundary: string;
  };
  fast_principal_cube: {
    affine_searches: number;
    completed_27_bit_cubes: number;
    completed_28_bit_cubes: number;
    completed_31_bit_cubes: number;
    completed_32_bit_cubes: number;
    exact_assignment_visits: string;
    assignments_per_engine_second: string;
    frontier_tie_assignments: number;
    proven_new_dephased_states_relative_prior: string;
    proven_new_dephased_states_across_pinned_audits: string;
    zero_pivot_corrections: number;
    new_local_h_classes: number;
    strict_wins: number;
    claim_boundary: string;
  };
  campaign_extensions: CampaignExtension[];
  diversified_search: {
    eo_exact_directions: string;
    gomea_exact_evaluations: string;
    gomea_archive_entries: number;
    tempering_exact_proposals: string;
    tempering_elite_reseeds: number;
    quality_diversity_exact_directions: string;
    quality_diversity_archive_cells: number;
    quality_diversity_frontier_cores_seen: number;
    new_h_only_classes: number;
    local_h_classes_after: number;
    local_ht_classes_after: number;
    strict_wins: number;
    claim_boundary: string;
  };
  gram_campaign_snapshot: {
    campaigns: number;
    exact_screens: string;
    stored_square_grams: string;
    hasse: {
      rejected: string;
      no_obstruction: string;
    };
    shell: {
      tested: string;
      rejected: string;
    };
    sign_factor_survivors: number;
    shell_milp: {
      signed_column_canonical_factors: number;
      frontier_ties: number;
      strict_wins: number;
    };
    claim_boundary: string;
  };
};

type CampaignSnapshotProps = {
  snapshot: Snapshot;
};

function formatInteger(value: string | number): string {
  return BigInt(value).toLocaleString("en-US");
}

function formatMetricValue(value: string): string {
  return /^\d+$/.test(value) ? formatInteger(value) : value;
}

function formatDate(value: string): string {
  return new Intl.DateTimeFormat("en-US", {
    day: "numeric",
    month: "short",
    year: "numeric",
    timeZone: "UTC",
  }).format(new Date(`${value}T00:00:00Z`));
}

export function CampaignSnapshot({ snapshot }: CampaignSnapshotProps) {
  const audits = snapshot.exact_audits;
  const portal = snapshot.portal_campaign;
  const neutral = snapshot.neutral_frontier;
  const qubo = snapshot.qubo_trust_region;
  const cubes = snapshot.fast_principal_cube;
  const diverse = snapshot.diversified_search;
  const gram = snapshot.gram_campaign_snapshot;
  const newestExtension = snapshot.campaign_extensions[0];

  return (
    <section className="campaign-section" aria-labelledby="campaign-title">
      <div className="wrap">
        <div className="campaign-heading">
          <div>
            <p className="eyebrow">
              Campaign snapshot · {formatDate(snapshot.updated)}
            </p>
            <h2 id="campaign-title">New matrices. Same frontier.</h2>
          </div>
          <p>
            <strong>No strict score improvement.</strong>{" "}
            {newestExtension ? (
              <>
                The newest extension, {newestExtension.kicker}, logged{" "}
                {formatMetricValue(newestExtension.headline.value)}{" "}
                {newestExtension.headline.label}. {newestExtension.summary}
              </>
            ) : (
              <>No campaign extension is recorded yet.</>
            )}
          </p>
        </div>

        <section
          className="portal-campaign"
          aria-labelledby="portal-campaign-title"
        >
          <div className="portal-campaign-copy">
            <p className="portal-campaign-kicker">{portal.method}</p>
            <h3 id="portal-campaign-title">Two new portals. Same frontier.</h3>
            <p>
              The substantive run found{" "}
              <strong>
                {formatInteger(portal.substantive_frontier_ties)} exact ties
              </strong>
              , but no determinant above the target.
            </p>
            <dl className="portal-run-metrics">
              <div>
                <dt>New local classes</dt>
                <dd>+{portal.new_local_h_ht_classes}</dd>
              </div>
              <div>
                <dt>Frontier ties</dt>
                <dd>{formatInteger(portal.substantive_frontier_ties)}</dd>
              </div>
              <div>
                <dt>Strict wins</dt>
                <dd>{portal.strict_wins}</dd>
              </div>
            </dl>
          </div>

          <div className="portal-campaign-visual">
            <div
              className="portal-atlas-transition"
              role="img"
              aria-label={`Frozen local H/HT frontier atlas grew from ${portal.local_h_ht_classes_before} to ${portal.local_h_ht_classes_after} classes`}
            >
              <div className="portal-atlas-state">
                <span>Frozen local atlas</span>
                <strong>{portal.local_h_ht_classes_before}</strong>
                <div className="portal-atlas-nodes" aria-hidden="true">
                  {Array.from(
                    { length: portal.local_h_ht_classes_before },
                    (_, index) => (
                      <i key={`before-${index}`} />
                    ),
                  )}
                </div>
              </div>
              <div className="portal-atlas-route" aria-hidden="true">
                <span>+{portal.new_local_h_ht_classes}</span>
              </div>
              <div className="portal-atlas-state portal-atlas-state-after">
                <span>Local atlas now</span>
                <strong>{portal.local_h_ht_classes_after}</strong>
                <div className="portal-atlas-nodes" aria-hidden="true">
                  {Array.from(
                    { length: portal.local_h_ht_classes_after },
                    (_, index) => (
                      <i
                        className={
                          index >= portal.local_h_ht_classes_before
                            ? "portal-atlas-node-new"
                            : undefined
                        }
                        key={`after-${index}`}
                      />
                    ),
                  )}
                </div>
              </div>
            </div>

            <div
              className="portal-radius-audits"
              aria-label="Exact radius-three portal closures"
            >
              {portal.exact_radius_3_audits.map((audit) => (
                <article key={audit.id}>
                  <div className="portal-radius-glyph" aria-hidden="true">
                    <i />
                  </div>
                  <div>
                    <span>{audit.label} · exact radius ≤ 3</span>
                    <strong>{formatInteger(audit.assignments)}</strong>
                    <small>
                      assignments · {audit.frontier_ties} ties ·{" "}
                      {audit.strict_wins} wins · {audit.status}
                    </small>
                  </div>
                </article>
              ))}
            </div>

            <div className="portal-connector-summary">
              <span>Exact portal connectors</span>
              <strong>
                {portal.exact_connector_cubes.count} × 2<sup>32</sup>
              </strong>
              <small>
                {formatInteger(portal.exact_connector_cubes.assignments)}{" "}
                assignments · {portal.exact_connector_cubes.frontier_ties}{" "}
                endpoint ties · {portal.exact_connector_cubes.strict_wins} wins
              </small>
              <em>{portal.exact_connector_cubes.result}</em>
            </div>

            <p className="portal-claim-boundary">{portal.claim_boundary}</p>
          </div>
        </section>

        <div className="campaign-audits" aria-label="Completed exact audits">
          {audits.map((audit) => (
            <article key={audit.label}>
              <span>{audit.label}</span>
              <strong>{formatInteger(audit.assignments)}</strong>
              <p>{audit.detail}</p>
            </article>
          ))}
        </div>

        <div
          className="campaign-extensions"
          aria-label="Latest completed search routes"
        >
          {snapshot.campaign_extensions.map((extension) => (
            <article key={extension.id}>
              <header>
                <span>{extension.kicker}</span>
                <div>
                  <strong>{formatMetricValue(extension.headline.value)}</strong>
                  <small>{extension.headline.label}</small>
                </div>
              </header>
              <dl>
                {extension.metrics.map((metric) => (
                  <div key={metric.label}>
                    <dt>{metric.label}</dt>
                    <dd>{formatMetricValue(metric.value)}</dd>
                  </div>
                ))}
              </dl>
              <p>{extension.summary}</p>
            </article>
          ))}
        </div>

        <div className="campaign-route-strip" aria-label="Completed campaigns">
          <article>
            <span>Exact affine cubes</span>
            <strong>{formatInteger(cubes.exact_assignment_visits)}</strong>
            <p>
              {cubes.completed_27_bit_cubes} × 2^27 +{" "}
              {cubes.completed_28_bit_cubes} × 2^28 +{" "}
              {cubes.completed_31_bit_cubes} × 2^31 +{" "}
              {cubes.completed_32_bit_cubes} × 2^32 visits ·{" "}
              {cubes.strict_wins} wins
            </p>
            <small>
              {formatInteger(
                cubes.proven_new_dephased_states_across_pinned_audits,
              )}{" "}
              distinct dephased states certified across pinned affine audits.{" "}
              {cubes.claim_boundary}
            </small>
          </article>
          <article>
            <span>Diversity archive</span>
            <strong>
              {formatInteger(diverse.quality_diversity_archive_cells)} niches
            </strong>
            <p>
              {formatInteger(diverse.quality_diversity_exact_directions)}{" "}
              directions · {diverse.strict_wins} wins
            </p>
            <small>{diverse.claim_boundary}</small>
          </article>
          <article className="neutral-result">
            <span>Neutral network</span>
            <strong>{neutral.new_raw_matrices} new matrices</strong>
            <p>
              {neutral.exact_cycles} exact {neutral.states_per_cycle}-cycles ·{" "}
              {neutral.neutral_edges} edges · Hamming{" "}
              {neutral.max_hamming_from_baseline}
            </p>
            <small>{neutral.claim_boundary}</small>
          </article>
          <article>
            <span>QUBO trust region</span>
            <strong>+{qubo.new_local_h_classes} H-class</strong>
            <p>
              {formatInteger(qubo.exact_block_evaluations)} exact blocks ·{" "}
              {qubo.frontier_ties} tie · {qubo.strict_wins} wins
            </p>
            <small>
              {qubo.local_h_classes_after} local classes total.{" "}
              {qubo.claim_boundary}
            </small>
          </article>
        </div>

        <div className="gram-snapshot">
          <div className="gram-snapshot-copy">
            <p className="eyebrow">
              Gram route · {gram.campaigns} campaigns
            </p>
            <h3>
              {formatInteger(gram.exact_screens)} exact screens. Zero
              survivors.
            </h3>
            <p>{gram.claim_boundary}</p>
          </div>

          <div className="gram-snapshot-grid">
            <article>
              <span>Exact screens</span>
              <strong>{formatInteger(gram.exact_screens)}</strong>
              <p>{gram.campaigns} campaigns</p>
            </article>
            <article>
              <span>Square Grams</span>
              <strong>{formatInteger(gram.stored_square_grams)}</strong>
              <p>stored</p>
            </article>
            <article>
              <span>Hasse gate</span>
              <strong>{formatInteger(gram.hasse.rejected)} rejected</strong>
              <p>{formatInteger(gram.hasse.no_obstruction)} no obstruction</p>
            </article>
            <article>
              <span>Shell gate</span>
              <strong>
                {formatInteger(gram.shell.rejected)}/
                {formatInteger(gram.shell.tested)}
              </strong>
              <p>rejected</p>
            </article>
            <article className="gram-outcome">
              <span>Sign factors</span>
              <strong>{gram.sign_factor_survivors}</strong>
              <p>survivors</p>
            </article>
            <article>
              <span>Shell MILP</span>
              <strong>
                {gram.shell_milp.signed_column_canonical_factors}
              </strong>
              <p>
                frontier ties · {gram.shell_milp.strict_wins} wins
              </p>
            </article>
          </div>
        </div>
      </div>
    </section>
  );
}
