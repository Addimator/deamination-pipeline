use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::env;


fn main() -> io::Result<()> {



    let args: Vec<String> = env::args().collect();
    // Überprüfe, ob genügend Argumente übergeben wurden
    if args.len() < 4 {
        eprintln!("Usage: {} <bedGraph_file_path> <bases_file_path> <output_file_path>", args[0]);
        return Err(io::Error::new(io::ErrorKind::Other, "Not enough arguments"));
    }

    // Die Argumente sind nullbasiert, args[0] ist der Name des Programms selbst
    let bedGraph_file_path = &args[1];
    let bases_file_path = &args[2];
    let output_file_path = &args[3];


    let mut meth_pos_forward: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut meth_pos_reverse: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut unmeth_pos_forward: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut unmeth_pos_reverse: HashMap<String, HashMap<String, usize>> = HashMap::new();

    let bedgraph_file = File::open(bedGraph_file_path).expect("Unable to open bedGraph file");
    let bases_file = File::open(bases_file_path).expect("Unable to open bases file");

    let bedgraph_reader = BufReader::new(bedgraph_file);
    let mut bases_reader = BufReader::new(bases_file).lines();

    for bed_line in bedgraph_reader.lines() {
        let bed_line = bed_line.expect("Error reading bedGraph line");
        let bed_fields: Vec<&str> = bed_line.split('\t').collect();

        // Extract relevant information from bedgraph
        let chrom = bed_fields[0].to_string();
        let methylation = bed_fields[3].parse::<usize>().expect("Invalid methylation value");

        // Get the two lines from bases.txt for the current position
        let bases_line_forward = bases_reader.next().expect("Error reading bases forward line").unwrap();
        let bases_line_reverse = bases_reader.next().expect("Error reading bases reverse line").unwrap();

        let bases_fields_forward: Vec<&str> = bases_line_forward.split('\t').collect();
        let bases_fields_reverse: Vec<&str> = bases_line_reverse.split('\t').collect();

        // Update the appropriate HashMaps based on methylation status and direction
        if methylation > 0 {
            update_base_counts(&mut meth_pos_forward, chrom.clone(), bases_fields_forward);
            update_base_counts(&mut meth_pos_reverse, chrom.clone(), bases_fields_reverse);
        } else {
            update_base_counts(&mut unmeth_pos_forward, chrom.clone(), bases_fields_forward);
            update_base_counts(&mut unmeth_pos_reverse, chrom.clone(), bases_fields_reverse);
        }
    }

    // Print the results
    println!("meth_pos_forward: {:?}", meth_pos_forward);
    println!("meth_pos_reverse: {:?}", meth_pos_reverse);
    println!("unmeth_pos_forward: {:?}", unmeth_pos_forward);
    println!("unmeth_pos_reverse: {:?}", unmeth_pos_reverse);
    Ok(())
}

fn update_base_counts(
    target_map: &mut HashMap<String, HashMap<String, usize>>,
    chrom: String,
    bases_fields: Vec<&str>,
) {
    let dir = bases_fields[2].to_string();
    let entry = target_map
        .entry(chrom.clone())
        .or_insert_with(|| HashMap::new());
    for (base, count_str) in bases_fields.iter().skip(3).enumerate() {
        let count = count_str.parse::<usize>().expect("Invalid base count");
        entry.insert(base.to_string(), count);
    }
}
