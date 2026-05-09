export default function TumorPhPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">Tumor pH Dashboard</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="heat" aria-label="reaction diffusion pH cross section" /></div>
      <section className="grid">
        <div className="panel"><h2>Warburg</h2><div className="metric">pH_e 7.4 to 6.7</div><div className="metric">CAIX inhibition raises pH</div></div>
        <div className="panel"><h2>CEST</h2><div className="metric">Sparse Z-spectrum to pH map</div><div className="metric">Tumor overlay</div></div>
      </section>
    </main>
  );
}
