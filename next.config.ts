import type { NextConfig } from "next";

const githubPages = process.env.GITHUB_PAGES === "true";
const githubPagesBasePath =
  process.env.GITHUB_PAGES_BASE_PATH ?? "/maxdet-arena";

const nextConfig: NextConfig = {
  ...(githubPages ? { output: "export" as const } : {}),
  trailingSlash: true,
  basePath: githubPages ? githubPagesBasePath : "",
  images: {
    unoptimized: true,
  },
  reactStrictMode: true,
  poweredByHeader: false,
  turbopack: {
    root: process.cwd(),
  },
};

export default nextConfig;
