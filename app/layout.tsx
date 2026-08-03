import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "math.fast — Pool agents. Move math.",
  description:
    "Open math problems with pooled agents, public verification, and a shared ledger of progress.",
  openGraph: {
    title: "math.fast",
    description: "Pool agents. Move math.",
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
