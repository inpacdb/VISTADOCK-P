import os
import multiprocessing
import argparse
import csv
from pipeline_worker_functions import run_screening_track, ensure_directory
from ligand_prep import run_ligand_prep
import sys

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()

# Mapping user selection to OpenMM Force FIelds
FORCE_FIELD_MAP = {
    "amber14-all": ["amber14-all.xml", "implicit/gbn2.xml"], #Default FF
    "amber14-SB": ["amber14/protein.ff14SB.xml", "amber14/DNA.OL15.xml", "amber14/RNA.OL3.xml", "implicit/gbn2.xml"],
    "amber19-SB": ["amber19/protein.ff19SB.xml", "implicit/gbn2.xml"],
    "amber99sb": ["amber99sb.xml", "implicit/obc2.xml"],
    "amber99sbildn": ["amber99sbildn", "implicit/obc2.xml"],
    "amber03": ["amber03.xml", "implicit/obc2.xml"],
    "amber10": ["amber10.xml", "implicit/obc2.xml"],
    "amber10": ["amber10.xml", "implicit/obc2.xml"],
    "charmm36": ["charmm36.xml", "implicit/obc2.xml"],
    "amoeba2013": ["amoeba2013.xml", "amoeba2013_gk.xml"],
    "amoeba2018": ["amoeba2018.xml", "amoeba2018_gk.xml"]
}

# FORMAT GBSA
def format_gbsa(val):
    """Helper to safely format GBSA values."""
    if isinstance(val, (int, float)):
        return f"{val:.2f}"
    return str(val)

def get_safe_name(name):
    return "".join([c for c in name if c.isalnum() or c in ('_','-')])

# GENERATE REPORTS 
def generate_reports(results_vina, results_vinardo, output_base, scoring_mode):
    """
    Generates three reports:
    1. A Text Table (.txt) for quick reading.
    2. A CSV File (.csv) for Excel/Data Analysis (Links OPEN the file).
    3. An HTML File (.html) for a Web-Like experience (Links DOWNLOAD the file).
    """
    all_ligands = sorted(set(results_vina.keys()) | set(results_vinardo.keys()))
    
    txt_file = os.path.join(output_base, "FINAL_SUMMARY.txt")
    csv_file = os.path.join(output_base, "FINAL_SUMMARY.csv")
    html_file = os.path.join(output_base, "FINAL_SUMMARY.html")
    
    # HTML HEADER & CSS
    html_content = ["""
    <!DOCTYPE html>
    <html>
    <head>
        <title>VISTADOCK-P Results</title>
        <style>
            body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; background-color: #f9f9f9; }
            h2 { color: #333; margin: 0; }
            table { border-collapse: collapse; width: 100%; box-shadow: 0 2px 8px rgba(0,0,0,0.1); background: white; }
            
            /* Header styling for sorting indicators */
            th { 
                background-color: #007bff; color: white; padding: 12px; text-align: left; 
                position: sticky; top: 0; cursor: pointer; user-select: none; padding-right: 25px;
            }
            th:hover { background-color: #0056b3; }
            th::after { content: '\\2195'; position: absolute; right: 8px; color: rgba(255,255,255,0.5); font-size: 0.9em; }
            th.sort-asc::after { content: '\\25B2'; color: white; }
            th.sort-desc::after { content: '\\25BC'; color: white; }
            
            td { border-bottom: 1px solid #ddd; padding: 8px; color: #333; }
            tr:hover { background-color: #f1f1f1; }
            .btn { 
                background-color: #28a745; color: white; padding: 6px 12px; 
                text-decoration: none; border-radius: 4px; font-size: 13px; font-weight: bold;
                display: inline-block;
            }
            .btn:hover { background-color: #218838; }
            .na { color: #ccc; font-style: italic; }
            
            /* NEW: Header layout and Export Button styling */
            .header-container { display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px; }
            .export-btn {
                background-color: #6f42c1; color: white; padding: 8px 16px; border: none;
                border-radius: 4px; font-size: 14px; font-weight: bold; cursor: pointer;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2);
            }
            .export-btn:hover { background-color: #59339d; }
        </style>
    </head>
    <body>
    <div class="header-container">
        <h2>VISTADOCK-P Final Summary</h2>
        <button id="exportHtmlBtn" class="export-btn">Save Current View as HTML</button>
    </div>
    <table>
    <thead>
        <tr>
            <th>Ligand</th>
            <th>Vina</th>
            <th>Vinardo</th>
            <th>CNN (V)</th>
            <th>CNN (VO)</th>
            <th>GBSA (V)</th>
            <th>GBSA (VO)</th>
            <th>PLIP (Vina)</th>
            <th>PLIP (Vinardo)</th>
            <th>Vina Complex</th>
            <th>Vinardo Complex</th>
        </tr>
    </thead>
    <tbody>
    """]

    # generate content
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # CSV Header
        header = [
            'LIGAND', 
            'VINA_SCORE', 'VINARDO_SCORE', 'SCORE_DIFF', 
            'CNN_SCORE_VINA', 'CNN_AFFINITY_VINA',
            'CNN_SCORE_VINARDO', 'CNN_AFFINITY_VINARDO',
            'GBSA_VINA', 'GBSA_VINARDO', 'GBSA_DIFF',
            'PLIP_VINA', 'PLIP_VINARDO', 
            'CONSENSUS_HIT',
            'LINK_VINA_COMPLEX', 'LINK_VINARDO_COMPLEX'
        ]
        writer.writerow(header)
        
        for lig in all_ligands:
            v_data = results_vina.get(lig, {'smina': "Excluded during ranking", 'gbsa': "N/A", 'plip': '-', 'properties': {}})
            vo_data = results_vinardo.get(lig, {'smina': "Excluded during ranking", 'gbsa': "N/A", 'plip': '-', 'properties': {}})
            
            s_v = v_data['smina']
            s_vo = vo_data['smina']
            
            if isinstance(s_v, (int, float)) and isinstance(s_vo, (int, float)):
                diff_smina = f"{abs(s_v - s_vo):.2f}"
            else:
                diff_smina = "N/A"
                
            str_s_v = f"{s_v:.2f}" if isinstance(s_v, (int, float)) else str(s_v)
            str_s_vo = f"{s_vo:.2f}" if isinstance(s_vo, (int, float)) else str(s_vo)
            
            g_v = v_data['gbsa']
            g_vo = vo_data['gbsa']
            if isinstance(g_v, (int, float)) and isinstance(g_vo, (int, float)):
                diff_gbsa = abs(g_v - g_vo)
            else:
                diff_gbsa = "N/A"
            
            is_consensus = (s_v != "Excluded during ranking" and s_vo != "Excluded during ranking")

            cnn_score_v = v_data.get('properties', {}).get('CNNscore', 'N/A')
            cnn_aff_v = v_data.get('properties', {}).get('CNNaffinity', 'N/A')
            cnn_score_vo = vo_data.get('properties', {}).get('CNNscore', 'N/A')
            cnn_aff_vo = vo_data.get('properties', {}).get('CNNaffinity', 'N/A')

            safe_name = get_safe_name(lig)

            if s_v != "Excluded during ranking":
                path_v = f"./vina_track/complexes/{safe_name}_complex.pdb" 
                link_v_csv = f'=HYPERLINK("{path_v}", "Open PDB")'
                link_v_html = f'<a href="{path_v}" class="btn" download="{safe_name}_vina.pdb">Download</a>'
            else:
                link_v_csv = "N/A"
                link_v_html = '<span class="na">N/A</span>'

            if s_vo != "Excluded during ranking":
                path_vo = f"./vinardo_track/complexes/{safe_name}_complex.pdb"
                link_vo_csv = f'=HYPERLINK("{path_vo}", "Open PDB")'
                link_vo_html = f'<a href="{path_vo}" class="btn" download="{safe_name}_vinardo.pdb">Download</a>'
            else:
                link_vo_csv = "N/A"
                link_vo_html = '<span class="na">N/A</span>'
            
            writer.writerow([
                lig, 
                str_s_v, str_s_vo, diff_smina,
                cnn_score_v, cnn_aff_v,
                cnn_score_vo, cnn_aff_vo,
                format_gbsa(g_v), format_gbsa(g_vo), format_gbsa(diff_gbsa),
                v_data['plip'], vo_data['plip'],
                str(is_consensus),
                link_v_csv, link_vo_csv
            ])

            html_content.append(f"""
            <tr>
                <td><b>{lig}</b></td>
                <td>{str_s_v}</td>
                <td>{str_s_vo}</td>
                <td>{cnn_score_v}</td>
                <td>{cnn_score_vo}</td>
                <td>{format_gbsa(g_v)}</td>
                <td>{format_gbsa(g_vo)}</td>
                <td style="font-size:0.9em">{v_data['plip']}
                <td style="font-size:0.9em; color:#666">{vo_data['plip']}</td>
                <td>{link_v_html}</td>
                <td>{link_vo_html}</td>
            </tr>
            """)
            
    html_content.append("""
    </tbody>
    </table>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            const table = document.querySelector('table');
            const headers = table.querySelectorAll('th');
            const tbody = table.querySelector('tbody');
            let currentSortCol = -1;
            let isAsc = true;

            headers.forEach((th, index) => {
                th.addEventListener('click', () => {
                    if (currentSortCol === index) {
                        isAsc = !isAsc;
                    } else {
                        currentSortCol = index;
                        isAsc = true;
                    }

                    headers.forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
                    th.classList.add(isAsc ? 'sort-asc' : 'sort-desc');

                    const rows = Array.from(tbody.querySelectorAll('tr'));
                    rows.sort((rowA, rowB) => {
                        let cellA = rowA.cells[index].innerText.trim();
                        let cellB = rowB.cells[index].innerText.trim();

                        const parseVal = (val) => {
                            if (val === 'N/A' || val === '-') return isAsc ? Infinity : -Infinity;
                            let num = parseFloat(val);
                            return isNaN(num) ? val : num;
                        };

                        let valA = parseVal(cellA);
                        let valB = parseVal(cellB);

                        if (typeof valA === 'number' && typeof valB === 'number') {
                            return isAsc ? valA - valB : valB - valA;
                        }
                        return isAsc ? cellA.localeCompare(cellB) : cellB.localeCompare(cellA);
                    });
                    tbody.append(...rows);
                });
            });

            const exportBtn = document.getElementById('exportHtmlBtn');
            if(exportBtn) {
                exportBtn.addEventListener('click', function() {
                    // Capture the entire live HTML document
                    let liveHtml = '<!DOCTYPE html>\\n' + document.documentElement.outerHTML;
                    
                    // Create a blob containing the text
                    const blob = new Blob([liveHtml], { type: 'text/html' });
                    const url = URL.createObjectURL(blob);
                    
                    // Create a temporary hidden link to trigger the download
                    const a = document.createElement('a');
                    a.href = url;
                    // Provide a default file name
                    a.download = 'VISTADOCK-P_Sorted_Summary.html'; 
                    document.body.appendChild(a);
                    a.click();
                    
                    // Clean up
                    document.body.removeChild(a);
                    URL.revokeObjectURL(url);
                });
            }
        });
    </script>
    </body>
    </html>
    """)
    
    with open(html_file, 'w') as f:
        f.write("\n".join(html_content))

    print(f"[REPORT] Reports generated:")
    print(f"         1. {csv_file}")
    print(f"         2. {html_file} (Use this for 'Download PDB' buttons)")

    # GENERATE TEXT TABLE
    w_lig = len("LIGAND")
    w_plip_v = len("PLIP (VINA)")
    w_plip_vo = len("PLIP (VINARDO)")

    for lig in all_ligands:
        w_lig = max(w_lig, len(lig))
        p_v = str(results_vina.get(lig, {}).get('plip', '-'))
        w_plip_v = max(w_plip_v, len(p_v))
        p_vo = str(results_vinardo.get(lig, {}).get('plip', '-'))
        w_plip_vo = max(w_plip_vo, len(p_vo))
    
    w_lig += 2; w_plip_v += 2; w_plip_vo += 2

    header_str = (
        f"{'LIGAND':<{w_lig}} | "
        f"{'VINA':<7} | {'VINARDO':<7} | {'DIFF(S)':<7} | "
        f"{'CNN (VINA)':<12} | {'CNN (VINARDO)':<12} | "
        f"{'GBSA(V)':<8} | {'GBSA(VO)':<8} | {'DIFF(G)':<7} | "
        f"{'PLIP (VINA)':<{w_plip_v}} | {'PLIP (VINARDO)':<{w_plip_vo}}"
    )
    divider = "-" * len(header_str)
    
    consensus_rows = []

    with open(txt_file, 'w') as f:
        f.write(header_str + "\n" + divider + "\n")
        print("\n" + header_str)
        print(divider)
        
        for lig in all_ligands:
            v_data = results_vina.get(lig, {'smina': "Excluded during ranking", 'gbsa': "N/A", 'plip': '-', 'properties': {}})
            vo_data = results_vinardo.get(lig, {'smina': "Excluded during ranking", 'gbsa': "N/A", 'plip': '-', 'properties': {}})
            
            s_v = v_data['smina']
            s_vo = vo_data['smina']
            
            if isinstance(s_v, (int, float)) and isinstance(s_vo, (int, float)):
                diff_smina = f"{-abs(s_v - s_vo):.2f}"
            else:
                diff_smina = "N/A"
                
            str_s_v = f"{s_v:.2f}" if isinstance(s_v, (int, float)) else str(s_v)
            str_s_vo = f"{s_vo:.2f}" if isinstance(s_vo, (int, float)) else str(s_vo)

            cnn_v = v_data['properties'].get('CNNscore', '-')
            if cnn_v != '-': cnn_v = f"{float(cnn_v):.2f}"
            cnn_vo = vo_data['properties'].get('CNNscore', '-')
            if cnn_vo != '-': cnn_vo = f"{float(cnn_vo):.2f}"

            g_v = v_data['gbsa']
            g_vo = vo_data['gbsa']
            
            if isinstance(g_v, (int, float)) and isinstance(g_vo, (int, float)):
                d_g = f"-{abs(g_v - g_vo):.2f}"
            else:
                d_g = "-"

            row = (
                f"{lig:<{w_lig}} | "
                f"{str_s_v:<7} | {str_s_vo:<7} | {diff_smina:<7} | "
                f"{cnn_v:<12} | {cnn_vo:<12} | "
                f"{format_gbsa(g_v):<8} | {format_gbsa(g_vo):<8} | {d_g:<7} | "
                f"{str(v_data['plip']):<{w_plip_v}} | {str(vo_data['plip']):<{w_plip_vo}}"
            )
            
            f.write(row + "\n")
            print(row)
            
            if s_v != "Excluded during ranking" and s_vo != "Excluded during ranking":
                consensus_rows.append(row)

        if consensus_rows:
            cons_header = "\n\n" + "="*40 + " CONSENSUS CANDIDATES " + "="*40
            f.write(cons_header + "\n")
            print(cons_header)
            f.write(header_str + "\n" + divider + "\n")
            print(header_str)
            print(divider)
            for row in consensus_rows:
                f.write(row + "\n")
                print(row)

if __name__ == '__main__':

    multiprocessing.set_start_method('spawn', force=True)
    pipeline_description = (
        "VISTADOCK-P: Virtual Integrative Screening & Triage Docking Pipeline\n"
        "-------------------------------------------------------------------------------------------------------------------------\n"
        "Developed by: Umashankar Vetrivel & Nidhish Kumar S.\n"
        "Bioinformatics Division, Indian Council of Medical Research - National Institute for Research in Tuberculosis (ICMR-NIRT)\n"
        "-------------------------------------------------------------------------------------------------------------------------"
    )
    
    parser = argparse.ArgumentParser(description=pipeline_description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    # Required Files
    parser.add_argument('--receptor_pdb', required=True, help='Path to receptor PDB file')
    parser.add_argument('--ligand', required=True, help='Path to ligand SDF file')
    parser.add_argument('--output_dir', default='pipeline_results', help='Directory to save results\n')

    # Box Coordinates
    parser.add_argument('--center_x', type=float, required=True, help='Center X')
    parser.add_argument('--center_y', type=float, required=True, help='Center Y')
    parser.add_argument('--center_z', type=float, required=True, help='Center Z\n')
    parser.add_argument('--size_x', type=float, required=True, help='Size X')
    parser.add_argument('--size_y', type=float, required=True, help='Size Y')
    parser.add_argument('--size_z', type=float, required=True, help='Size Z\n')

    # Lig prep params
    parser.add_argument('--skip_prep', action='store_true', help="Skip ligand preparation step")
    parser.add_argument('--prep_ff', type=str, default="MMFF94", help="Force Field for ligand preparation")

    parser.add_argument('--cpu_count', type=int, default=0, help='Total CPUs to allocate (Selecting 0 causes the tool to utilises max. no. of CPUs available)')
    parser.add_argument('--prep_cpus', type=int, default=0, help='CPUs for ligand preparation (0 = Auto/75%%)')

    # Workflow Params
    parser.add_argument('--exh_rapid', type=int, default=8, help='Exhaustiveness RAPID (def: 8)')
    parser.add_argument('--exh_balanced', type=int, default=16, help='Exhaustiveness BALANCED (def: 16)')
    parser.add_argument('--exh_ultra', type=int, default=32, help='Exhaustiveness ULTRA (def: 32)\n')
    
    parser.add_argument('--frac_rapid', type=float, default=50, help='Percentage of retained top ligands RAPID (def: 50%)')
    parser.add_argument('--frac_balanced', type=float, default=30, help='Percentage of retained top ligands BALANCED (def: 30%)')
    parser.add_argument('--frac_ultra', type=float, default=100, help='Percentage of retained top ligands ULTRA (def: 100%)\n')
    
    parser.add_argument('--scoring', choices=['vina', 'vinardo', 'both'], default='both', help='Scoring function (def: both)')
    parser.add_argument('--no_mmgbsa', action='store_true', help='Skip MM-GBSA')
    parser.add_argument('--single_step', action='store_true', help="Bypass RAPID & BALANCED steps, run only one docking step")
    parser.add_argument('--no_plip', action='store_true', help='Skip PLIP Analysis')
    parser.add_argument('--lipinski', action='store_true', help="Filter input ligands using Lipinski's Rule of 5")
    parser.add_argument('--cnn_scoring', type=str, default='rescore', help="Enable GNINA CNN Re-Scoring (Requires GPU)")
    parser.add_argument('--gpu_device', type=int, default=0, help="GPU Device ID to use (def: 0)")

    # MMGBSA params
    parser.add_argument('--forcefield', default='amber14-all', choices=FORCE_FIELD_MAP.keys(), help='Force Field for MM-GBSA')
    
    parser.add_argument('--gbsa_temp', type=float, default=300.0, help='Temperature in Kelvin for MM-GBSA (def: 300.0 K)')
    parser.add_argument('--gbsa_friction', type=float, default=1.0, help='Co-efficient of Friction in 1/picosecond for MM-GBSA (def: 1.0/ps)')
    parser.add_argument('--gbsa_timestep', type=float, default=2.0, help='Timestep in femtoseconds for MM-GBSA (def: 2.0 fs)')

    args = parser.parse_args()

    max_sys_cpus = multiprocessing.cpu_count()
    user_cpus = args.cpu_count

    if user_cpus > max_sys_cpus or user_cpus < 0:
        print(f"[CONFIGURATION] Requested CPUs ({user_cpus}) invalid or greater than the number of CPUs available ({max_sys_cpus}). Using maximum number of CPUs available.")
        total_cpus = max_sys_cpus
    else:
        total_cpus = user_cpus

    scoring_mode=args.scoring
    if scoring_mode == 'both':
        cpus_per_track = max(1, int(total_cpus//2))
        print(f"[CONFIGURATION] Scoring Mode: BOTH. Allocating {cpus_per_track} CPUs per track (Total: {cpus_per_track*2})")
    else:
        cpus_per_track = total_cpus
        print (f"[CONFIGURATION] Scoring Mode: {scoring_mode.upper()}. Allocating {cpus_per_track} CPUs to track.")

    # SETUP
    output_base = args.output_dir
    ensure_directory(output_base)
    
    log_file_path = os.path.join(output_base, "vistadockp.log")
    sys.stdout = Logger(log_file_path)
    sys.stderr = sys.stdout
    
    print("\n" + "="*50)
    print("VISTADOCK-P EXECUTION LOG STARTED")
    print("="*50 + "\n")
    
    config = vars(args)

    config['cpu_count'] = cpus_per_track

    selected_ff = args.forcefield
    config['ff-xmls'] = FORCE_FIELD_MAP.get(selected_ff, FORCE_FIELD_MAP['amber14-all'])
    print(f"[CONFIGURATION] Force Field: {selected_ff} : {config['ff-xmls']}")  

    # LIGAND PREP
    print ("===== LIGAND PREPARATION =====")
    
    if args.skip_prep:
        print("Skipping ligand preparation (User requested).")
        prepared_lig_path = args.ligand
    else:
        prepared_lig_path = run_ligand_prep(
            input_file=args.ligand,
            output_dir=args.output_dir,
            apply_lipinski=args.lipinski,
            cpu_count=args.prep_cpus,
            prep_ff = args.prep_ff
        )
    
    if not prepared_lig_path:
        print("[CRITICAL ERROR] Ligand preparation failed. Exiting....")
        exit(1)

    config['ligand'] = prepared_lig_path
    
    manager = multiprocessing.Manager()
    final_results = manager.dict()
    processes = []

    scoring_mode = args.scoring
    modes_to_run = ['vina', 'vinardo'] if scoring_mode == 'both' else [scoring_mode]
        
    print(f"===== LAUNCHING PIPELINE: {args.scoring.upper()} =====")
    print(f"Receptor: {args.receptor_pdb}")
    print(f"Ligand:   {args.ligand}")
    
    for mode in modes_to_run:
        p = multiprocessing.Process(
            target=run_screening_track, 
            args=(mode, output_base, final_results, config)
        )
        processes.append(p)
        p.start()
        
    for p in processes:
        p.join()
        
    print("\n===== ALL TRACKS COMPLETED. GENERATING REPORTS =====")
    
    # Report generation based on mode
    if scoring_mode == 'both':
        generate_reports(final_results['vina'], final_results['vinardo'], output_base, 'both')
    elif scoring_mode == 'vina':
        generate_reports(final_results['vina'], {}, output_base, 'vina')
    else:
        generate_reports({}, final_results['vinardo'], output_base, 'vinardo')
        
    print(f"[Done] Pipeline Finished.")