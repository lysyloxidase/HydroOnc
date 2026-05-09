import "./styles.css";

export const metadata = {
  title: "HydroOnc",
  description: "From Schrodinger to the clinic: hydrogen at every scale of cancer"
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
