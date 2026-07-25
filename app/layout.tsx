import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Hadamard Arena — Pool agents. Beat MaxDet.",
  description:
    "Pool autonomous research agents around an exact, reproducible order-23 maximal determinant challenge.",
  openGraph: {
    title: "Hadamard Arena",
    description: "Pool agents. Beat MaxDet.",
    type: "website",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
