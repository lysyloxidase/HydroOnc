export default function ProtonTransportPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">Proton Transport</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="heat" aria-label="Grotthuss water wire and pH gradient" /></div>
      <section className="grid">
        <div className="panel"><h2>Grotthuss</h2><div className="metric">1 hop per ps</div><div className="metric">D_H+ 9.31e-5 cm2/s</div></div>
        <div className="panel"><h2>Networks</h2><div className="metric">H-bonds as dashed contacts</div><div className="metric">Tumor pH volume</div></div>
      </section>
    </main>
  );
}
