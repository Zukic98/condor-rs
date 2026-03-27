// labrador/examples/poc.rs
use std::time::Instant;

use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::witness::Witness;
use labrador::relation::statement::Statement;
use labrador::commitments::common_instances::AjtaiInstances;
use labrador::prover::LabradorProver;
// ADDED: Import the verifier
use labrador::verifier::LabradorVerifier;
use labrador::transcript::sponges::shake::ShakeSponge;
use labrador::transcript::LabradorTranscript;

fn main() {
    println!("==================================================");
    println!("🚀 LaBRADOR Falcon Aggregation - Proof of Concept");
    println!("==================================================\n");

    println!("[*] Loading environment and parameters...");
    let ep = EnvironmentParameters::default();

    println!("[*] Generating Witness (Simulating Falcon signatures)...");
    let witness = Witness::new(ep.rank, ep.multiplicity, ep.beta_sq);

    println!("[*] Generating public Statement...");
    let st = Statement::new(&witness, &ep);

    println!("[*] Generating CRS matrices (AjtaiInstances)...");
    let crs = AjtaiInstances::new(&ep);

    println!("\n>>> STARTING PROVER (Aggregation Process) <<<\n");
    let start_prove = Instant::now();

    let mut prover = LabradorProver::new(&ep, &crs, &witness, &st);
    
    // The proof variable is now passed to the verifier, resolving previous warnings.
    let proof: LabradorTranscript<ShakeSponge> = prover.prove().expect("Error generating proof!");

    let prove_time = start_prove.elapsed();
    println!("[+] Aggregated proof generated in: {:?}", prove_time);

    println!("\n>>> STARTING VERIFIER (Proof Verification) <<<\n");
    let start_verify = Instant::now();

    // Initialize the verifier (it does not need the witness, as it has no access to secret data/signatures)
    let mut verifier = LabradorVerifier::new(&ep, &crs, &st);
    
    // Verifier checks the proof
    let verify_result = verifier.verify(&proof);
    
    let verify_time = start_verify.elapsed();
    
    match verify_result {
        Ok(_) => println!("[+] SUCCESS! Proof is VALID. (Verification time: {:?})", verify_time),
        Err(e) => println!("[-] ERROR! Proof is invalid: {:?}", e),
    }

    println!("\n==================================================");
    println!("🎉 Full LaBRADOR cycle (Aggregation + Verification) complete!");
    println!("==================================================");
}