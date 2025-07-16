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
    
    def create_reference_mapping(self):
        """
        Create mapping of which domain should be used as reference for each domain
        """
        connectivity = self.build_connectivity_graph()
        
        for i, domain_id in enumerate(self.domain_hierarchy):
            if i == 0:
                # First domain is global reference
                self.reference_mapping[domain_id] = domain_id
            else:
                # Look for direct connections to earlier domains in hierarchy
                # Search in REVERSE order of hierarchy to find most recently processed neighbor
                best_reference = None
                for j in range(i - 1, -1, -1):  # Go backwards through earlier domains
                    earlier_domain_id = self.domain_hierarchy[j]
                    if earlier_domain_id in connectivity[domain_id]:
                        best_reference = earlier_domain_id
                        break  # Take the first (most recent) connection found
                
                # Fallback to global reference if no direct connection
                if best_reference is None:
                    best_reference = self.domain_hierarchy[0]
                
                self.reference_mapping[domain_id] = best_reference
        
        return self.reference_mapping
    
    def get_analysis_pairs(self):
        """
        Get list of (moving_domain, reference_domain) pairs for analysis
        """
        hierarchy = self.create_hierarchical_ordering()
        reference_mapping = self.create_reference_mapping()
        
        analysis_pairs = []
        
        for domain_id in hierarchy:
            reference_id = reference_mapping[domain_id]
            
            if domain_id != reference_id:  # Skip self-references
                analysis_pairs.append((domain_id, reference_id))
        print(f"Analysis pairs: {analysis_pairs}")
        return analysis_pairs
    
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