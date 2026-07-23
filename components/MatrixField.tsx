type MatrixFieldProps = {
  matrix: number[][];
};

export function MatrixField({ matrix }: MatrixFieldProps) {
  return (
    <div
      className="matrix-frame"
      role="img"
      aria-label="Published order-23 sign matrix"
    >
      <div className="matrix-grid" aria-hidden="true">
        {matrix.flatMap((row, rowIndex) =>
          row.map((entry, columnIndex) => (
            <span
              className={entry === 1 ? "matrix-plus" : "matrix-minus"}
              key={`${rowIndex}-${columnIndex}`}
              title={`row ${rowIndex + 1}, column ${columnIndex + 1}: ${entry}`}
            />
          )),
        )}
      </div>
      <div className="matrix-axis matrix-axis-x" aria-hidden="true">
        23 columns
      </div>
      <div className="matrix-axis matrix-axis-y" aria-hidden="true">
        23 rows
      </div>
    </div>
  );
}
