import pyLCIO
import numpy as np
import uproot

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open("/local/d1/mu+mu-/sim/1000_10_sim.slcio")
for ievt, event in enumerate(reader):
    print("Event", ievt)
    mcpCollection = event.getCollection("MCParticle")
    for mcp in mcpCollection:
        vertex = mcp.getVertex()
        r_xy = np.sqrt(vertex[0]**2 + vertex[1]**2)
        z = vertex[2]
        id = mcp.id()
        if abs(mcp.getPDG())==15:
            print("\tTau Parent - unique id, PDG:", [(x.id(), x.getPDG()) for x in mcp.getParents()] if mcp.getParents() else None, "Tau vertex, status:", [r_xy,z], mcp.getGeneratorStatus())
        endp = mcp.getEndpoint()
        r_xy = np.sqrt(endp[0]**2 + endp[1]**2)
        z = endp[2]
        if abs(mcp.getPDG())==1000015 or abs(mcp.getPDG())==2000015:
            print("\tStau - unique id, status, endpoint:", id, mcp.getGeneratorStatus(), [r_xy,z])
                     
    
    for coll_name, (sim_x, sim_y, sim_z, sim_time, sim_pdg, sim_mcpid, sim_layer), (reco_x, reco_y, reco_z, reco_time, reco_layer) in [
                ("VertexBarrelCollection", ("sim_VB_x", "sim_VB_y", "sim_VB_z", "sim_VB_time", "sim_VB_pdg", "sim_VB_mcpid", "sim_VB_layer"), ("reco_VB_x", "reco_VB_y", "reco_VB_z", "reco_VB_time", "reco_VB_layer")),
                ("VertexEndcapCollection", ("sim_VE_x", "sim_VE_y", "sim_VE_z", "sim_VE_time", "sim_VE_pdg", "sim_VE_mcpid", "sim_VE_layer"), ("reco_VE_x", "reco_VE_y", "reco_VE_z", "reco_VE_time", "reco_VE_layer")),
                ("InnerTrackerBarrelCollection", ("sim_IB_x", "sim_IB_y", "sim_IB_z", "sim_IB_time", "sim_IB_pdg", "sim_IB_mcpid", "sim_IB_layer"), ("reco_IB_x", "reco_IB_y", "reco_IB_z", "reco_IB_time", "reco_IB_layer")),
                ("InnerTrackerEndcapCollection", ("sim_IE_x", "sim_IE_y", "sim_IE_z", "sim_IE_time", "sim_IE_pdg", "sim_IE_mcpid", "sim_IE_layer"), ("reco_IE_x", "reco_IE_y", "reco_IE_z", "reco_IE_time", "reco_IE_layer")),
                ("OuterTrackerBarrelCollection", ("sim_OB_x", "sim_OB_y", "sim_OB_z", "sim_OB_time", "sim_OB_pdg", "sim_OB_mcpid", "sim_OB_layer"), ("reco_OB_x", "reco_OB_y", "reco_OB_z", "reco_OB_time", "reco_OB_layer")),
                ("OuterTrackerEndcapCollection", ("sim_OE_x", "sim_OE_y", "sim_OE_z", "sim_OE_time", "sim_OE_pdg", "sim_OE_mcpid", "sim_OE_layer"), ("reco_OE_x", "reco_OE_y", "reco_OE_z", "reco_OE_time", "reco_OE_layer"))
        ]:       
        collection = event.getCollection(coll_name)
        for hit in collection:
                # print("hit",collection.getTypeName(), collection.getParameters())
                position = hit.getPosition()
                mcp = hit.getMCParticle()
                hit_pdg = mcp.getPDG() if mcp else np.nan
                if np.abs(hit_pdg) == 1000015 or np.abs(hit_pdg) == 2000015:
                    position = hit.getPosition()
                    r_xy = np.sqrt(position[0]**2 + position[1]**2)
                    z = position[2]
                    print("\tStau hit detected", [r_xy,z])
