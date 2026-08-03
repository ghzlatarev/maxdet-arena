import type { Metadata } from "next";

export const metadata: Metadata = {
  title: "Solved problems — math.fast",
  description:
    "A curated public record of research problems solved with substantial agent involvement.",
};

export default function SolvedLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return children;
}
