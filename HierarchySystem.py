import numpy as np
from copy import deepcopy

class HierarchicalDomainSystem:
    """
    Implements a Fortran DynDom-style hierarchical domain reference system
    where domains are processed in connectivity order and each serves as
    a local reference point for its analysis.
    """
    
    def __init__(self, domains):
        self.domains = domains
        self.domain_hierarchy = []
        self.reference_mapping = {}
        
    def build_connectivity_graph(self):
        """
        Build connectivity graph showing which domains are directly connected
        """
        connectivity = {}
        
        for domain in self.domains:
            connectivity[domain.domain_id] = []
            
        # Check connectivity between all domain pairs
        for i, domain_a in enumerate(self.domains):
            for j, domain_b in enumerate(self.domains):
                if i == j:
                    continue
                    
                if self.are_domains_connected(domain_a, domain_b):
                    connectivity[domain_a.domain_id].append(domain_b.domain_id)
                    
        return connectivity
    
    def are_domains_connected(self, domain_a, domain_b):
        """
        Check if two domains are directly connected (adjacent in sequence)
        """
        # Get all segment boundaries for both domains
        for seg_a in domain_a.segments:
            for seg_b in domain_b.segments:
                # Check if domains are adjacent (one segment ends where another begins)
                if seg_a[1] + 1 == seg_b[0] or seg_b[1] + 1 == seg_a[0]:
                    return True
        return False
    
    def create_hierarchical_ordering(self):
        """
        Create hierarchical domain ordering following Fortran DynDom logic
        FIXED: Better tie-breaking when domains have equal connectivity
        """
        connectivity = self.build_connectivity_graph()
        remaining_domains = set(domain.domain_id for domain in self.domains)
        
        print("Initial connectivity:", connectivity)
        
        hierarchy = []
        
        while remaining_domains:
            # Find domains with maximum connectivity among remaining domains
            max_connectivity = -1
            candidates = []
            
            for domain_id in remaining_domains:
                current_connections = [conn for conn in connectivity[domain_id] 
                                    if conn in remaining_domains]
                conn_count = len(current_connections)
                
                if conn_count > max_connectivity:
                    max_connectivity = conn_count
                    candidates = [domain_id]
                elif conn_count == max_connectivity:
                    candidates.append(domain_id)
            
            # TIE-BREAKING: If multiple domains have same connectivity
            if len(candidates) > 1:
                print(f"Tie between domains {candidates} with connectivity {max_connectivity}")
                
                # Tie-breaker: Choose domain with most residues
                domain_sizes = {}
                for domain_id in candidates:
                    domain = self.domains[domain_id]
                    domain_sizes[domain_id] = sum(domain.segments[:, 1] + 1 - domain.segments[:, 0])
                
                most_connected_domain = max(candidates, key=lambda d: domain_sizes[d])
                
                print(f"Tie-breaker: Selected domain {most_connected_domain} (size: {domain_sizes[most_connected_domain]} residues)")
                
            else:
                most_connected_domain = candidates[0]
            
            print(f"Selected domain {most_connected_domain} with connectivity {max_connectivity}")
            
            # Add to hierarchy
            hierarchy.append(most_connected_domain)
            
            # Remove this domain from remaining set
            remaining_domains.remove(most_connected_domain)
            
            # Remove this domain from all other domains' connectivity lists
            for domain_id in remaining_domains:
                if most_connected_domain in connectivity[domain_id]:
                    connectivity[domain_id].remove(most_connected_domain)
            
            print(f"Updated connectivity after removing {most_connected_domain}:", 
                {k: v for k, v in connectivity.items() if k in remaining_domains})
        
        self.domain_hierarchy = hierarchy
        print(f"Final domain hierarchy: {self.domain_hierarchy}")
        return hierarchy
    
    def get_analysis_pairs(self):
        """
        Get list of (moving_domain, reference_domain) pairs for analysis
        """
        # MUST create hierarchy first, then reference mapping
        hierarchy = self.create_hierarchical_ordering()
        reference_mapping = self.create_reference_mapping()
        
        analysis_pairs = []
        
        for domain_id, reference_list in reference_mapping.items():
            if isinstance(reference_list, list):
                for reference_id in reference_list:
                    if domain_id != reference_id:  # Skip self-references
                        analysis_pairs.append((domain_id, reference_id))
            else:
                # Handle case where reference_mapping still returns single values
                reference_id = reference_list
                if domain_id != reference_id:
                    analysis_pairs.append((domain_id, reference_id))
        
        print(f"Analysis pairs: {analysis_pairs}")
        return analysis_pairs

    def create_reference_mapping(self):
        """
        Create mapping - small domains analyzed relative to larger connected domains
        BUT respect the global reference constraint
        """
        # Make sure hierarchy exists
        if not self.domain_hierarchy:
            self.create_hierarchical_ordering()
        
        connectivity = self.build_connectivity_graph()
        self.reference_mapping = {}
        
        # Get global reference from hierarchy
        global_reference_id = self.domain_hierarchy[0]
        
        # Get domain sizes
        domain_sizes = {}
        for domain in self.domains:
            domain_sizes[domain.domain_id] = domain.num_residues
        
        for domain in self.domains:
            domain_id = domain.domain_id
            
            if domain_id == global_reference_id:
                # Global reference is never a moving domain
                self.reference_mapping[domain_id] = [domain_id]  # Self-reference
                continue
            
            # Find all connected domains that are larger OR are the global reference
            reference_candidates = []
            for connected_id in connectivity[domain_id]:
                if (domain_sizes[connected_id] > domain_sizes[domain_id] or 
                    connected_id == global_reference_id):
                    reference_candidates.append(connected_id)
            
            # ONLY analyze if there are valid reference candidates
            if reference_candidates:
                self.reference_mapping[domain_id] = reference_candidates
            else:
                # Don't create any analysis pairs for this domain
                self.reference_mapping[domain_id] = [domain_id]  # Self-reference (no analysis)
        
        return self.reference_mapping
    
    def print_analysis_plan(self):
        """
        Print the analysis plan showing domain hierarchy and reference mapping
        """
        hierarchy = self.create_hierarchical_ordering()
        reference_mapping = self.create_reference_mapping()
        
        print("\n=== HIERARCHICAL DOMAIN ANALYSIS PLAN ===")
        print("Domain Hierarchy (processing order):", hierarchy)
        print("\nReference Mapping:")
        
        for domain_id in hierarchy:
            reference_id = reference_mapping[domain_id]
            if domain_id == reference_id:
                print(f"  Domain {domain_id}: GLOBAL REFERENCE (stationary)")
            else:
                print(f"  Domain {domain_id}: analyzed relative to Domain {reference_id}")
        
        print("\nAnalysis pairs:")
        analysis_pairs = self.get_analysis_pairs()
        for moving, reference in analysis_pairs:
            print(f"  Domain {moving} movement → relative to Domain {reference}")