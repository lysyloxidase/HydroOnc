export default function HBondPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">Oncoprotein H-Bonds</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="heat" aria-label="protein H-bond difference map" /></div>
      <section className="grid">
        <div className="panel"><h2>p53</h2><div className="metric">R175H loses 4+ H-bonds</div><div className="metric">R273H DNA contact below 10%</div></div>
        <div className="panel"><h2>GNN</h2><div className="metric">PDBbind test affinity</div><div className="caveat">APR-246 Phase 3 failed.</div></div>
      </section>
    </main>
  );
}
