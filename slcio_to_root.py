import pyLCIO
import ROOT
from ROOT import std, TFile, TTree
import glob
from typing import Dict
import numpy as np
from math import *
import multiprocessing
import os
import re

############################################################################################################
'''
Run the script with:
nohup python slcio_to_root.py > reco_v2/output.log 2>&1 &
'''

# Set the maximum number of events to process
max_events = -1 # Set to -1 to run over all events
output_directory = "/local/d1/lrozanov/mucoll-tutorial-2023/reco_v2"

# Gather input files, excluding certain patterns if necessary
fnames = glob.glob("/local/d1/mu+mu-/reco_v2/*.slcio") 
exclude_patterns = ['_light', 'ns']
# Filter filenames: include only those that do not contain any of the exclusion patterns
fnames = [fname for fname in fnames if not any(excl in fname for excl in exclude_patterns)]

print("Found %i files." % len(fnames))

def initialize_vector_and_branch(tree, branches: Dict[str, str], nested_branch_names: set):
    """
    Creates vectors and branches in a tree for given names and data types.

    Parameters:
    - tree: The ROOT TTree object to which the branches will be added.
    - branches: A dictionary where keys are the branch names and values are the data types ('float' or 'int').
    """
    vectors = {}
    for name, dtype in branches.items():
        if dtype not in ['float', 'int']:
            raise ValueError("Unsupported data type: {}".format(dtype))
        if name in nested_branch_names:
            vec = std.vector("std::vector<{}>".format(dtype))()
        else:
            vec = std.vector(dtype)()
        tree.Branch(name, vec)
        vectors[name] = vec
    return vectors

def clear_vectors(vectors, nested_branch_names):
    """
    Clears all vectors stored in the provided dictionary, with specific handling for nested vectors.

    Parameters:
    - vectors: A dictionary of vectors to be cleared.
    - nested_branch_names: A set of branch names that are nested and need special clearing.
    """
    for name, vector in vectors.items():
        if name in nested_branch_names:
            # Clear each inner vector in the nested vector
            for v in vector:
                v.clear()
            vector.clear()
        else:
            vector.clear()

def strip_filename(filename, output_directory):
    pattern = '.slcio'
    base_filename = os.path.basename(filename)
    # Use re.sub to remove '_sim.slcio' or '_digi.slcio' from the end of the filename
    root_filename = re.sub(pattern, '.root', base_filename)
    output_path = os.path.join(output_directory, root_filename)
    return output_path

# Define the branches and their data types
branches_info = {
    # MCPs
    "mcp_pt": 'float', "mcp_phi": 'float', "mcp_eta": 'float', "pdgid": 'int', "status": 'int', "prod_vertex_x": 'float', 
    "prod_vertex_y": 'float', "prod_vertex_z": 'float', "prod_time": 'float',  "id": 'int', 
    # Staus
    "mcp_stau_pt": 'float', "mcp_stau_eta": 'float', "mcp_stau_phi": 'float',
    # Tracks
    "d0": 'float', "z0": 'float', "nhits": 'int', "pixel_nhits": 'int', "inner_nhits": 'int', "outer_nhits": 'int',
    "track_pt": 'float', "track_eta": 'float', "track_phi": 'float', "ndf": 'int', "chi2": 'float',
    # Track hits (outer)
    "x": 'float', "y": 'float', "z": 'float', "time": 'float', "corrected_time": 'float', 
    "hit_layer": 'int', "hit_detector": 'int', "hit_side": 'int',
    # Sim hits
    "sim_VB_x": 'float', "sim_VB_y": 'float', "sim_VB_z": 'float', "sim_VB_time": 'float', "sim_VB_pdg": 'float', "sim_VB_mcpid": 'float', "sim_VB_layer": 'float',
    "sim_VE_x": 'float', "sim_VE_y": 'float', "sim_VE_z": 'float', "sim_VE_time": 'float', "sim_VE_pdg": 'float', "sim_VE_mcpid": 'float', "sim_VE_layer": 'float',
    "sim_IB_x": 'float', "sim_IB_y": 'float', "sim_IB_z": 'float', "sim_IB_time": 'float', "sim_IB_pdg": 'float', "sim_IB_mcpid": 'float', "sim_IB_layer": 'float',
    "sim_IE_x": 'float', "sim_IE_y": 'float', "sim_IE_z": 'float', "sim_IE_time": 'float', "sim_IE_pdg": 'float', "sim_IE_mcpid": 'float', "sim_IE_layer": 'float',
    "sim_OB_x": 'float', "sim_OB_y": 'float', "sim_OB_z": 'float', "sim_OB_time": 'float', "sim_OB_pdg": 'float', "sim_OB_mcpid": 'float', "sim_OB_layer": 'float',
    "sim_OE_x": 'float', "sim_OE_y": 'float', "sim_OE_z": 'float', "sim_OE_time": 'float', "sim_OE_pdg": 'float', "sim_OE_mcpid": 'float', "sim_OE_layer": 'float',
    # Reco hits
    "reco_VB_x": 'float', "reco_VB_y": 'float', "reco_VB_z": 'float', "reco_VB_time": 'float', "reco_VB_layer": 'float',
    "reco_VE_x": 'float', "reco_VE_y": 'float', "reco_VE_z": 'float', "reco_VE_time": 'float', "reco_VE_layer": 'float',
    "reco_IB_x": 'float', "reco_IB_y": 'float', "reco_IB_z": 'float', "reco_IB_time": 'float', "reco_IB_layer": 'float',
    "reco_IE_x": 'float', "reco_IE_y": 'float', "reco_IE_z": 'float', "reco_IE_time": 'float', "reco_IE_layer": 'float',
    "reco_OB_x": 'float', "reco_OB_y": 'float', "reco_OB_z": 'float', "reco_OB_time": 'float', "reco_OB_layer": 'float',
    "reco_OE_x": 'float', "reco_OE_y": 'float', "reco_OE_z": 'float', "reco_OE_time": 'float', "reco_OE_layer": 'float'    
}

# Handle hits as nested branches
nested_branch_names = {"x", "y", "z", "time", "corrected_time", "hit_layer", "hit_detector", "hit_side"}

speedoflight = 299792458/1000000  # mm/ns

def process_file(slcio_file):
    # Set up some options
    num_events = max_events # Set to -1 to run over all events
    n_mcp_stau = 0
    stau_hits = 0
    
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(slcio_file)
    
    # Create a new ROOT file and tree
    stripped_filename = strip_filename(slcio_file, output_directory)
    root_file = TFile(stripped_filename, "RECREATE")
    tree = TTree("Events", "LLPs Data")

    # Initialize vectors and branches
    vectors = initialize_vector_and_branch(tree, branches_info, nested_branch_names)   

    # Event loop starts here
    for ievt, event in enumerate(reader):
        if num_events > 0 and ievt >= num_events: break
        if ievt % 10 == 0: print("Processing event %i." % ievt)
        clear_vectors(vectors, nested_branch_names)
        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        trackCollection = event.getCollection("SiTracks")

        hit_collections = []
        IBTrackerHits = event.getCollection('ITBarrelHits')
        hit_collections.append(IBTrackerHits)
        IETrackerHits = event.getCollection('ITEndcapHits')
        hit_collections.append(IETrackerHits)
        OBTrackerHits = event.getCollection('OTBarrelHits')
        hit_collections.append(OBTrackerHits)
        OETrackerHits = event.getCollection('OTEndcapHits')
        hit_collections.append(OETrackerHits)
        VBTrackerHits = event.getCollection('VXDBarrelHits')
        hit_collections.append(VBTrackerHits)
        VETrackerHits = event.getCollection('VXDEndcapHits')
        hit_collections.append(VETrackerHits)
        
        # Relations
        VBrelationCollection = event.getCollection('VXDBarrelHitsRelations')
        VBrelation = pyLCIO.UTIL.LCRelationNavigator(VBrelationCollection)
        
        VErelationCollection = event.getCollection('VXDEndcapHitsRelations')
        VErelation = pyLCIO.UTIL.LCRelationNavigator(VErelationCollection)
        
        IBrelationCollection = event.getCollection('ITBarrelHitsRelations')
        IBrelation = pyLCIO.UTIL.LCRelationNavigator(IBrelationCollection)
        
        IErelationCollection = event.getCollection('ITEndcapHitsRelations')
        IErelation = pyLCIO.UTIL.LCRelationNavigator(IErelationCollection)
        
        OBrelationCollection = event.getCollection('OTBarrelHitsRelations')
        OBrelation = pyLCIO.UTIL.LCRelationNavigator(OBrelationCollection)
        
        OErelationCollection = event.getCollection('OTEndcapHitsRelations')
        OErelation = pyLCIO.UTIL.LCRelationNavigator(OErelationCollection)
        
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            
            vectors["mcp_pt"].push_back(mcp_tlv.Perp())
            vectors["mcp_phi"].push_back(mcp_tlv.Phi())
            vectors["mcp_eta"].push_back(mcp_tlv.Eta())
            vectors["pdgid"].push_back(mcp.getPDG())
            vectors["status"].push_back(mcp.getGeneratorStatus())
            vectors["prod_vertex_x"].push_back([mcp.getVertex()[i] for i in range(3)][0])
            vectors["prod_vertex_y"].push_back([mcp.getVertex()[i] for i in range(3)][1])
            vectors["prod_vertex_z"].push_back([mcp.getVertex()[i] for i in range(3)][2])
            vectors["prod_time"].push_back(mcp.getTime())
            vectors["id"].push_back(mcp.id())
            
            # Find the stau MCPs
            if abs(mcp.getPDG())==1000015 or abs(mcp.getPDG())==2000015:
                vectors["mcp_stau_pt"].push_back(mcp_tlv.Perp())
                vectors["mcp_stau_phi"].push_back(mcp_tlv.Phi())
                vectors["mcp_stau_eta"].push_back(mcp_tlv.Eta())
                n_mcp_stau += 1
        
        for track in trackCollection:
            Bfield = 5 #T, 3.57 for legacy
            theta = np.pi/2- np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            pt  = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            track_tlv = ROOT.TLorentzVector()
            track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
            
            vectors["d0"].push_back(track.getD0())
            vectors["z0"].push_back(track.getZ0())
            vectors["nhits"].push_back(track.getTrackerHits().size())
            vectors["track_pt"].push_back(pt)
            vectors["track_eta"].push_back(eta)
            vectors["track_phi"].push_back(phi)
            vectors["ndf"].push_back(track.getNdf())
            vectors["chi2"].push_back(track.getChi2())
            
            # Create local vectors for track hits
            local_vectors = {
                "track_x": std.vector('float')(),
                "track_y": std.vector('float')(),
                "track_z": std.vector('float')(),
                "track_time": std.vector('float')(),
                "track_corrected_time": std.vector('float')(),
                "track_hit_layer": std.vector('int')(),
                "track_hit_detector": std.vector('int')(),
                "track_hit_side": std.vector('int')()
            }
            
            pixel_nhit = 0
            inner_nhit = 0
            outer_nhit = 0
            for hit in track.getTrackerHits():
                # now decode hits
                encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                decoder = pyLCIO.UTIL.BitField64(encoding)
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                detector = decoder["system"].value()
                layer = decoder['layer'].value()
                side = decoder["side"].value()
                if detector == 1 or detector == 2:
                     pixel_nhit += 1
                if detector == 3 or detector == 4:
                    inner_nhit += 1
                if detector == 5 or detector == 6:
                    outer_nhit += 1
                    
                position = hit.getPosition()
                
                d = sqrt(position[0]*position[0] + position[1]
                     * position[1] + position[2]*position[2])
                tof = d/speedoflight

                resolution = 0.03
                if detector > 2:
                    resolution = 0.06

                corrected_time = hit.getTime()*(1.+ROOT.TRandom3(ievt).Gaus(0., resolution)) - tof
                
                 # Append hit data to local vectors
                local_vectors["track_x"].push_back(position[0])
                local_vectors["track_y"].push_back(position[1])
                local_vectors["track_z"].push_back(position[2])
                local_vectors["track_time"].push_back(hit.getTime())
                local_vectors["track_corrected_time"].push_back(corrected_time)
                local_vectors["track_hit_layer"].push_back(layer)
                local_vectors["track_hit_detector"].push_back(detector)
                local_vectors["track_hit_side"].push_back(side)

            # Append each track's vectors to the corresponding nested vector in the main dictionary
            for key in local_vectors:
                vectors[key.replace('track_', '')].push_back(local_vectors[key])
            
        for coll_name, relation, (sim_x, sim_y, sim_z, sim_time, sim_pdg, sim_mcpid, sim_layer), (reco_x, reco_y, reco_z, reco_time, reco_layer) in [
            ("VertexBarrelCollection", VBrelation, ("sim_VB_x", "sim_VB_y", "sim_VB_z", "sim_VB_time", "sim_VB_pdg", "sim_VB_mcpid", "sim_VB_layer"), ("reco_VB_x", "reco_VB_y", "reco_VB_z", "reco_VB_time", "reco_VB_layer")),
            ("VertexEndcapCollection", VErelation, ("sim_VE_x", "sim_VE_y", "sim_VE_z", "sim_VE_time", "sim_VE_pdg", "sim_VE_mcpid", "sim_VE_layer"), ("reco_VE_x", "reco_VE_y", "reco_VE_z", "reco_VE_time", "reco_VE_layer")),
            ("InnerTrackerBarrelCollection", IBrelation, ("sim_IB_x", "sim_IB_y", "sim_IB_z", "sim_IB_time", "sim_IB_pdg", "sim_IB_mcpid", "sim_IB_layer"), ("reco_IB_x", "reco_IB_y", "reco_IB_z", "reco_IB_time", "reco_IB_layer")),
            ("InnerTrackerEndcapCollection", IErelation, ("sim_IE_x", "sim_IE_y", "sim_IE_z", "sim_IE_time", "sim_IE_pdg", "sim_IE_mcpid", "sim_IE_layer"), ("reco_IE_x", "reco_IE_y", "reco_IE_z", "reco_IE_time", "reco_IE_layer")),
            ("OuterTrackerBarrelCollection", OBrelation, ("sim_OB_x", "sim_OB_y", "sim_OB_z", "sim_OB_time", "sim_OB_pdg", "sim_OB_mcpid", "sim_OB_layer"), ("reco_OB_x", "reco_OB_y", "reco_OB_z", "reco_OB_time", "reco_OB_layer")),
            ("OuterTrackerEndcapCollection", OErelation, ("sim_OE_x", "sim_OE_y", "sim_OE_z", "sim_OE_time", "sim_OE_pdg", "sim_OE_mcpid", "sim_OE_layer"), ("reco_OE_x", "reco_OE_y", "reco_OE_z", "reco_OE_time", "reco_OE_layer"))
        ]:       
            try:
                collection = event.getCollection(coll_name)
                # THIS SEEMS TO SAVE BIB HITS TOO
                # can try some comparison between sim mcpid and mcpid to not save bib hits
                for hit in collection:
                    # print("hit",collection.getTypeName(), collection.getParameters())
                    position = hit.getPosition()
                    vectors[sim_x].push_back(position[0])
                    vectors[sim_y].push_back(position[1])
                    vectors[sim_z].push_back(position[2])
                    vectors[sim_time].push_back(hit.getTime())
                    
                    mcp = hit.getMCParticle()
                    hit_pdg = mcp.getPDG() if mcp else np.nan
                    mcpid = mcp.id() if mcp else np.nan
                    vectors[sim_pdg].push_back(hit_pdg)
                    vectors[sim_mcpid].push_back(mcpid)
                    
                    # Get layer
                    encoding = collection.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    detector = decoder["system"].value()
                    layer = decoder['layer'].value()
                    side = decoder["side"].value()
                    vectors[sim_layer].push_back(layer)

                    # Get the corresponding reco hit
                    reco_hit = relation.getRelatedFromObjects(hit)
                    if len(reco_hit) > 0:
                        recohit = reco_hit[0]
                        rel_position = recohit.getPosition()
                        
                        vectors[reco_x].push_back(rel_position[0])
                        vectors[reco_y].push_back(rel_position[1])
                        vectors[reco_z].push_back(rel_position[2])
                        vectors[reco_time].push_back(recohit.getTime())
                        vectors[reco_layer].push_back(layer)
                        
                        if np.abs(hit_pdg) == 1000015 or np.abs(hit_pdg) == 2000015:
                            stau_hits += 1
                    else:
                        vectors[reco_x].push_back(np.nan)
                        vectors[reco_y].push_back(np.nan)
                        vectors[reco_z].push_back(np.nan)
                        vectors[reco_time].push_back(np.nan)
                        vectors[reco_layer].push_back(np.nan)
       
            except Exception as e:
                print(f"Error accessing {coll_name}: {e}")
                        
        tree.Fill()
        ievt += 1
        
    root_file.Write()
    root_file.Close()
    reader.close()
    
    print(f"\nSummary statistics for {os.path.basename(slcio_file)}:")
    print("Ran over %i events."%ievt)
    print("Found:")
    print("\t%i Stau MCPs"%n_mcp_stau)
    print("\t%i Stau hits"%stau_hits)


n_cpu = len(fnames) if len(fnames) <= 16 else 16
pool = multiprocessing.Pool(n_cpu)

# Process each input file
pool.map(process_file, fnames)
# Close the pool of processes
pool.close()
pool.join()