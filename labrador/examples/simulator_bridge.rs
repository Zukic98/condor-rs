// labrador/examples/simulator_bridge.rs
use std::time::Instant;

use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::witness::Witness;
use labrador::relation::statement::Statement;
use labrador::commitments::common_instances::AjtaiInstances;
use labrador::prover::LabradorProver;
use labrador::verifier::LabradorVerifier;
use labrador::transcript::sponges::shake::ShakeSponge;
use labrador::transcript::LabradorTranscript;

fn main() {
    // 1. Initialize environment parameters and mock data
    let ep = EnvironmentParameters::default();
    let witness = Witness::new(ep.rank, ep.multiplicity, ep.beta_sq);
    let st = Statement::new(&witness, &ep);
    let crs = AjtaiInstances::new(&ep);

    // 2. Execute and measure the Prover (Aggregation)
    let start_prove = Instant::now();
    let mut prover = LabradorProver::new(&ep, &crs, &witness, &st);
    let proof: LabradorTranscript<ShakeSponge> = prover.prove().expect("Proof generation failed");
    
    // Convert duration to floating-point milliseconds for precise Python plotting
    let prove_time_ms = start_prove.elapsed().as_secs_f64() * 1000.0;

    // 3. Execute and measure the Verifier
    let start_verify = Instant::now();
    let mut verifier = LabradorVerifier::new(&ep, &crs, &st);
    let verify_result = verifier.verify(&proof);
    
    let verify_time_ms = start_verify.elapsed().as_secs_f64() * 1000.0;
    let is_valid = verify_result.is_ok();

    // 4. Output ONLY valid JSON so the Python simulator can parse it easily
    println!(
        r#"{{"prove_time_ms": {:.4}, "verify_time_ms": {:.4}, "is_valid": {}}}"#,
        prove_time_ms, verify_time_ms, is_valid
    );
}