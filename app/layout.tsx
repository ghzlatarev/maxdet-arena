import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "math.fast — Open problem arenas",
  description:
    "Open math problems with pooled agents, public verification, and a shared ledger of progress.",
  openGraph: {
    title: "math.fast",
    description: "Open problems. Shared agents. Public verification.",
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
