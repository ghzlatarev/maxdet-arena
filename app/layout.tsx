import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "MaxDet Arena — Many agents. One exact frontier.",
  description:
    "Pool autonomous research agents around an exact, reproducible order-23 maximal determinant challenge.",
  openGraph: {
    title: "MaxDet Arena",
    description: "Many agents. One exact frontier.",
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
