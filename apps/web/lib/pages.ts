export const pages = [
  { slug: "quantum", title: "Quantum Explorer", checks: ["1s spherical", "2p dumbbell", "3d multi-lobed"] },
  { slug: "proton-transport", title: "Proton Transport", checks: ["Grotthuss hopping", "pH gradient"] },
  { slug: "tumor-ph", title: "Tumor pH Dashboard", checks: ["Warburg pH_e 6.7", "CEST overlay"] },
  { slug: "proton-therapy", title: "Proton Therapy", checks: ["150 MeV at 15.6 cm", "SOBP flat"] },
  { slug: "hbond", title: "Oncoprotein H-Bonds", checks: ["p53 R175H lost bonds", "KRAS His95"] },
  { slug: "h2-therapy", title: "H2 Therapy", checks: ["GSH caveat", "anti-PD-1 synergy"] }
];

export const caveats = [
  "H2 selective scavenging is contested.",
  "APR-246 Phase 3 failed.",
  "Variable RBE is not ICRU standard.",
  "FLASH FAST-01 was feasibility, not efficacy.",
  "H2 cancer trials are small and not Phase 3.",
  "Research platform, not clinical decision support."
];
