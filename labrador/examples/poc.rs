// labrador/examples/poc.rs
use std::time::Instant;

use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::witness::Witness;
use labrador::relation::statement::Statement;
// Ovdje je ispravljena putanja koju nam je Rust sugerirao!
use labrador::commitments::common_instances::AjtaiInstances;
use labrador::prover::LabradorProver;
use labrador::transcript::sponges::shake::ShakeSponge;
use labrador::transcript::LabradorTranscript;

fn main() {
    println!("==================================================");
    println!("🚀 LaBRADOR Falcon Agregacija - Proof of Concept");
    println!("==================================================\n");

    println!("[*] Učitavanje okruženja i parametara (EnvironmentParameters)...");
    let ep = EnvironmentParameters::default();

    println!("[*] Generiranje 'Witness-a' (Simulacija Falcon potpisa)...");
    let witness = Witness::new(ep.rank, ep.multiplicity, ep.beta_sq);

    println!("[*] Generiranje javne tvrdnje (Statement)...");
    let st = Statement::new(&witness, &ep);

    println!("[*] Generiranje CRS matrica (AjtaiInstances)...");
    let crs = AjtaiInstances::new(&ep);

    println!("\n>>> POKREĆEM PROVER (Proces agregacije) <<<\n");
    
    let start_prove = Instant::now();
    
    let mut prover = LabradorProver::new(&ep, &crs, &witness, &st);
    let proof: LabradorTranscript<ShakeSponge> = prover.prove().expect("Greška prilikom generiranja dokaza!");
    
    let prove_time = start_prove.elapsed();
    
    println!("[+] USPJEH! Agregirani dokaz je generiran.");
    println!("[+] Vrijeme generiranja agregiranog dokaza: {:?}", prove_time);
    println!("[+] Veličina dokaza je kriptografski sažeta.\n");

    println!("==================================================");
    println!("Završeno! Sljedeći korak je dodavanje Verifikatora.");
    println!("==================================================");
}
