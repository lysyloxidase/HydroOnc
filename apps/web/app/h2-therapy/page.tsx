export default function H2TherapyPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">H2 Therapy</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="heat" aria-label="ROS and immune ODE curves" /></div>
      <section className="grid">
        <div className="panel"><h2>Immune ODE</h2><div className="metric">H2 + anti-PD-1 beats monotherapy</div><div className="metric">T_eff recovery</div></div>
        <div className="panel"><h2>Caveats</h2><div className="caveat">H2 selective scavenging is contested.</div><div className="caveat">H2 cancer trials are small and not Phase 3.</div></div>
      </section>
    </main>
  );
}
