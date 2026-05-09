export default function ProtonTherapyPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">Proton Therapy</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="heat" aria-label="Bragg peak and RBE heatmap" /></div>
      <section className="grid">
        <div className="panel"><h2>Bragg</h2><div className="metric">150 MeV peak at 15.6 cm</div><div className="metric">SOBP flat within 5%</div></div>
        <div className="panel"><h2>Caveats</h2><div className="caveat">Variable RBE is not ICRU standard.</div><div className="caveat">FLASH FAST-01 was feasibility, not efficacy.</div></div>
      </section>
    </main>
  );
}
