use std::ffi::c_float;
use std::fs::File;
use std::io::{BufReader, Result, BufRead, Write, BufWriter};
use structopt::StructOpt;
use std::path::PathBuf;
use find_bases::{extract_vcf_positions, count_bases_in_reads};
use assign_bases::update_base_counts;
use std::collections::HashMap;


mod find_bases;
mod assign_bases;


#[derive(Debug, StructOpt, Clone)]
pub enum Deamination {
    #[structopt(
        name = "find-bases",
        about = "Find bases at CpG sites",
        usage = "cargo run -- filter-bases aligned_reads.sam candidates.vcf pos_to_bases.txt"
    )]
    BaseFinder {
        #[structopt(
            name = "aligned-reads",
            parse(from_os_str),
            required = true,
            help = "Aligned reads in SAM format"
        )]
        sam_file_path: PathBuf,
        #[structopt(
            name = "candidates",
            parse(from_os_str),
            required = true,
            help = "Candidates of CpG positions in BCF format"
        )]
        bcf_file_path: PathBuf,
        #[structopt(
            name = "output", 
            parse(from_os_str), 
            help = "Path of output TXT file, if not given output is printed to stdout")]
        output: Option<PathBuf>,
    },
    #[structopt(
        name = "assign-bases",
        about = "Assign bases to methylated and unmethylated cases",
        usage = "cargo run -- assign-bases ref.bedGraph pos_to_bases.txt output.txt"
    )]
    BaseAssigner {
        #[structopt(
            name = "ref-bedgraph",
            parse(from_os_str),
            required = true,
            help = "Bedgraph for reference which positions are (un)methylated"
        )]
        bedGraph_path: PathBuf,
        #[structopt(
            name = "pos-to-bases",
            parse(from_os_str),
            required = true,
            help = "Bases on each CpG position"
        )]
        bases_file_path: PathBuf,
        #[structopt(
            name = "output", 
            parse(from_os_str), 
            help = "Path of output TXT file, if not given output is printed to stdout")]
        output: Option<PathBuf>,
    },
}




fn main() -> Result<()> {
    let opt = Deamination::from_args();
    match opt {
        Deamination::BaseFinder { sam_file_path, bcf_file_path, output } => {

            let vcf_positions = extract_vcf_positions(bcf_file_path)?;
            let position_counts = count_bases_in_reads(sam_file_path, &vcf_positions)?;

            // Öffnen Sie die Ausgabedatei zum Schreiben
            let output_file = File::create(output.unwrap())?;
            let mut writer = BufWriter::new(output_file);

            // Schreiben Sie die Ergebnisse in die Datei
            writeln!(&mut writer, "#CHROM	#POS	#DIR	#A	#C	#G	#T	#N")?;
            for ((reference_name, position, direction), counts) in &position_counts {
                write!(
                    &mut writer,
                    "{}	{}	{}	",
                    reference_name, position, direction
                )?;
                write!(&mut writer, "{}	{}	{}	{}	{}", counts.get(&'A').unwrap_or(&0), counts.get(&'C').unwrap_or(&0), counts.get(&'G').unwrap_or(&0), counts.get(&'T').unwrap_or(&0), counts.get(&'N').unwrap_or(&0))?;
                println!("{}	{}	{}	{}	{}", counts.get(&'A').unwrap_or(&0), counts.get(&'C').unwrap_or(&0), counts.get(&'G').unwrap_or(&0), counts.get(&'T').unwrap_or(&0), counts.get(&'N').unwrap_or(&0));
                writeln!(&mut writer)?;
            }
        }
        Deamination::BaseAssigner { bedGraph_path, bases_file_path, output } => {
            let mut meth_pos_forward: HashMap<String, HashMap<char, usize>> = HashMap::new();
            let mut meth_pos_reverse: HashMap<String, HashMap<char, usize>> = HashMap::new();
            let mut unmeth_pos_forward: HashMap<String, HashMap<char, usize>> = HashMap::new();
            let mut unmeth_pos_reverse: HashMap<String, HashMap<char, usize>> = HashMap::new();

            let bedgraph_file = File::open(bedGraph_path).expect("Unable to open bedGraph file");
            let bases_file = File::open(bases_file_path).expect("Unable to open bases file");

            let bedgraph_reader = BufReader::new(bedgraph_file);
            let bases_reader = BufReader::new(bases_file);
      




            let mut forward_baseline_map: HashMap<(String, usize), Vec<String>> = HashMap::new();
            let mut reverse_baseline_map: HashMap<(String, usize), Vec<String>> = HashMap::new();
            



            for base_line in bases_reader.lines() {
                let base_line = base_line.expect("Error reading bedGraph line");
                if base_line.starts_with("#") {
                    continue;
                }
                let base_fields: Vec<String>  = base_line.split('\t').map(|s| s.to_string()).collect();

                let chrom = base_fields[0].to_string();
                let position = base_fields[1].parse::<usize>().expect("Invalid position value");
                let direction = base_fields[2].parse::<char>().expect("Invalid direction");
            
                if direction == 'f'{
                    forward_baseline_map.insert((chrom.clone(), position), base_fields.clone());
                }
                else {
                    reverse_baseline_map.insert((chrom.clone(), position), base_fields.clone());
                }
            }

            for bed_line in bedgraph_reader.lines() {
                let bed_line = bed_line.expect("Error reading bedGraph line");
                if bed_line.starts_with("track") {
                    continue;
                }
                let bed_fields: Vec<&str> = bed_line.split('\t').collect();
            
                // Extrahiere relevante Informationen aus dem Bedgraph
                let chrom = bed_fields[0].to_string();
                let position = (bed_fields[1].parse::<usize>().expect("Invalid position value") + bed_fields[1].parse::<usize>().expect("Invalid position value")) / 2;
                let methylation = bed_fields[3].parse::<f64>().expect("Invalid methylation value");
                
                let chrom_clone1 = chrom.clone();
                let chrom_clone2 = chrom.clone();
                let chrom_clone3 = chrom.clone();
                let chrom_clone4 = chrom.clone();
                let chrom_clone5 = chrom.clone();
                let chrom_clone6 = chrom.clone();
                let chrom_clone7 = chrom.clone();
                if forward_baseline_map.contains_key(&(chrom, position)) {
                    let bases_fields_forward = &forward_baseline_map[&(chrom_clone1, position)];
                    let bases_fields_reverse = &reverse_baseline_map[&(chrom_clone2, position)];

                    if methylation > 20.0 {
                        update_base_counts(&mut meth_pos_forward, chrom_clone3, bases_fields_forward.to_vec());
                        update_base_counts(&mut meth_pos_reverse, chrom_clone4, bases_fields_reverse.to_vec());
                    } else {
                        update_base_counts(&mut unmeth_pos_forward, chrom_clone5, bases_fields_forward.to_vec());
                        update_base_counts(&mut unmeth_pos_reverse, chrom_clone6, bases_fields_reverse.to_vec());
                    }
                }
                else {
                    println!("Entry not found: ({}, {})", chrom_clone7, position)
                }

            }

            let output: Option<File> = match output {
                Some(path) => {
                    // If `output` contains a path, open the file for writing
                    Some(File::create(&path)?)
                }
                None => None, // If `output` is None, don't open any file
            };
        
            //Create a BCF writer depending on the output (to file or to stdout)
            let mut writer: Box<dyn Write> = match output {
                Some(file) => Box::new(BufWriter::new(file)),
                None => Box::new(BufWriter::new(std::io::stdout())),
            };

            // Print the results
            writeln!(writer, "meth_pos_forward: {:?}", meth_pos_forward);
            writeln!(writer, "meth_pos_reverse: {:?}", meth_pos_reverse);
            writeln!(writer, "unmeth_pos_forward: {:?}", unmeth_pos_forward);
            writeln!(writer, "unmeth_pos_reverse: {:?}", unmeth_pos_reverse);
            writeln!(writer, );

            for chrom in meth_pos_forward.keys(){
                writeln!(writer, "Methylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&meth_pos_forward, chrom, 'T'), bases_percentage(&meth_pos_reverse, chrom, 'A'));
                writeln!(writer, "Unmethylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&unmeth_pos_forward, chrom, 'T'), bases_percentage(&unmeth_pos_reverse, chrom, 'A'));
            }

        }
    }
    Ok(())
}

pub fn bases_percentage(target_map: &HashMap<String, HashMap<char, usize>>, chrom: &String, base: char) -> f64 {
    if let Some(map_chrom) = target_map.get(chrom) {
        let sum_bases: usize = map_chrom.values().sum();
        if let Some(&number_base) = map_chrom.get(&base) {
            if sum_bases > 0 {
                return (number_base as f64 / sum_bases as f64) * 100.0;
            }
        }
    }
    0.0 // Return 0.0 if something goes wrong or if sum_bases is 0.
}