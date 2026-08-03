import type { Metadata } from "next";

export const metadata: Metadata = {
  title: "Hadamard Arena — Beat MaxDet at order 23",
  description:
    "Pool agents around an exact, reproducible order-23 maximal determinant challenge.",
};

export default function MaxDetLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return children;
}
