export default function QuantumPage() {
  return (
    <main className="shell">
      <div className="topbar"><div className="brand">Quantum Explorer</div><a href="/">HydroOnc</a></div>
      <div className="viz"><div className="orbital" aria-label="1s 2p 3d orbital viewer" /></div>
      <section className="grid">
        <div className="panel"><h2>Orbitals</h2><div className="metric">1s spherical</div><div className="metric">2p dumbbell</div><div className="metric">3d multi-lobed</div></div>
        <div className="panel"><h2>Spectrum</h2><div className="metric">Balmer H-alpha 656.28 nm</div><div className="metric">Lyman alpha 121.57 nm</div></div>
      </section>
    </main>
  );
}
