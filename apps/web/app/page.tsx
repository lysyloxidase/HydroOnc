import { caveats, pages } from "../lib/pages";

export default function Home() {
  return (
    <main className="shell">
      <header className="topbar">
        <div>
          <div className="brand">HydroOnc</div>
          <div>From Schrodinger to the clinic: hydrogen at every scale of cancer</div>
        </div>
        <nav className="nav">
          {pages.map((page) => <a key={page.slug} href={`/${page.slug}`}>{page.title}</a>)}
        </nav>
      </header>
      <section className="grid">
        {pages.map((page) => (
          <article className="panel" key={page.slug}>
            <h2>{page.title}</h2>
            {page.checks.map((check) => <div className="metric" key={check}>{check}</div>)}
          </article>
        ))}
      </section>
      <section className="panel" style={{ marginTop: 18 }}>
        <h2>Caveats</h2>
        {caveats.map((caveat) => <div className="caveat" key={caveat}>{caveat}</div>)}
      </section>
    </main>
  );
}
