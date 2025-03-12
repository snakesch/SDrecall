import logging
import numpy as np

from src.log import logger

class HOMOSEQ_REGION:
    def __init__(self, vertex, graph, logger=logger):
        node_tuple = graph.vertex_properties["node_index"][vertex]
        self.chrom = node_tuple[0]
        self.start = node_tuple[1]
        self.end = node_tuple[2]
        self.strand = node_tuple[3]
        self.size = self.end - self.start
        # The subset corresponding with upstream region
        self.ups_rela_start = 0
        self.ups_rela_end = self.end - self.start
        # The subset corresponding with downstream region
        self.down_rela_start = 0
        self.down_rela_end = 0
        # The subset of the overlapping part between the upstream segment and the downstream corresponding segment
        self.rela_start = self.ups_rela_start
        self.rela_end = self.ups_rela_end
        self.data = node_tuple
        self.vertex = int(vertex) # here vertex should only be an integer ID, ensure the object is pickable for multiprocessing
        # rela_coords should be a dictionary that stores:
        # overlapping node (vertex) -> (rela_start, rela_end)
        # traverse_route looks like [(node, edge_type), (node, edge_type), (node, edge_type)]
        self.traverse_route = []

    def __hash__(self):
        return hash(self.data) + hash(self.vertex) + hash(self.ups_rela_start) + hash(self.ups_rela_end) + hash(self.down_rela_start) + hash(self.down_rela_end) + hash(self.rela_start) + hash(self.rela_end)

    def __repr__(self):
        return ", ".join(str(x) for x in [self.chrom, self.start, self.end, self.strand, self.size])

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end and self.strand == other.strand

    def __iter__(self):
        for item in [self.chrom, self.start + self.ups_rela_start, self.start + self.ups_rela_end, self.strand]:
            yield item

    def __getitem__(self, index):
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self.data))
            return [self.data[i] for i in range(start, stop, step)]
        else:
            return self.data[index]

    def __len__(self):
        return len(self.data)

    def fix_coord(self):
        return tuple([self.chrom, self.start + self.ups_rela_start, self.start + self.ups_rela_end, self.strand])
    
    def qnode_relative_region(self, qnode_tuple, logger=logger):
        # The qnode tuple should be a 4-item tuple containint (chrom, start, end, strand)
        # We also need to utilize the record info in the traverse route.
        total_route = list(self.traverse_route)
        assert qnode_tuple in [ n[0].data for n in total_route ], f"The qnode tuple {qnode_tuple} is not in the traverse route of the node {self.data}: {total_route}"
        qnode_route = total_route[total_route.index([ t for t in total_route if t[0].data == qnode_tuple ][0]):]
        if total_route != qnode_route:
            logger.info(f"The total route is \n{total_route}, and the qnode route is \n{qnode_route}\n")
        
        # Start to process the relative coordinates of the self node corresponding to the input qnode_tuple
        current_node = self
        qnode_rela_start = current_node.rela_start
        qnode_rela_end = current_node.rela_end
        assert qnode_rela_end > qnode_rela_start, f"The current node {current_node} does not have a valid rela_start {qnode_rela_start} and rela_end {qnode_rela_end}"

        while len(qnode_route) > 0:
            # Inside this loop, we calculate the relative start and end coordinates of the upstream node in the traverse route
            # Use 
            if type(qnode_route) == tuple:
                upstream_qnode, ups_edg_type = qnode_route
            else:
                upstream_qnode, ups_edg_type = qnode_route.pop() # Here upstream_qnode is a HOMOSEQ_REGION object
            if ups_edg_type == "segmental_duplication":
                if upstream_qnode.strand == current_node.strand:
                    qnode_rela_start = max(qnode_rela_start, 0)
                    qnode_rela_end = min(qnode_rela_end, upstream_qnode.end - upstream_qnode.start)
                else:
                    qsize = qnode_rela_end - qnode_rela_start
                    qnode_rela_end = min(current_node.end - current_node.start - qnode_rela_start, upstream_qnode.end - upstream_qnode.start)
                    qnode_rela_start = max(qnode_rela_end - qsize, 0)
            else:
                qnode_abs_start = qnode_rela_start + current_node.start
                qnode_abs_end = qnode_rela_end + current_node.start
                qsize = qnode_rela_end - qnode_rela_start
                
                current_rela_start = qnode_rela_start
                current_rela_end = qnode_rela_end

                qnode_rela_start = qnode_abs_start - upstream_qnode.start
                qnode_rela_end = min(qnode_rela_start + qsize, upstream_qnode.end - upstream_qnode.start)
                qnode_rela_start = max(qnode_rela_start, 0)

            if qnode_rela_end - qnode_rela_start <= 0:
                logger.warning(f"The relative end is smaller than the relative start for {current_node}, for current node, before updating the rela_start is {current_rela_start} and rela_end is {current_rela_end}. For the upstream qnode {upstream_qnode} (relationship: {ups_edg_type}), the rela_start calculated is {qnode_rela_start} and rela_end calculated at this iteration is {qnode_rela_end}. The original rela_start is {upstream_qnode.rela_start} and rela_end is {upstream_qnode.rela_end}")
                return "NaN", "NaN"

            current_node = upstream_qnode

        return qnode_rela_start, qnode_rela_end

# - Possible optimization - #   
# import logging
# import numpy as np

# logger = logging.getLogger("SDrecall")

# class HOMOSEQ_REGION:
#     def __init__(self, vertex, graph):
#         chrom, start, end, strand = graph.vertex_properties["node_index"][vertex]
#         self.chrom = chrom
#         self.start = start
#         self.end = end
#         self.strand = strand
#         self.size = end - start
#         self.ups_rela_start = 0
#         self.ups_rela_end = self.size
#         self.down_rela_start = 0  # Currently unused, consider removing
#         self.down_rela_end = 0    # Currently unused, consider removing
#         self.rela_start = 0  # Simplified: Initial overlap is with itself
#         self.rela_end = self.size # Simplified: Initial overlap is with itself
#         self.data = (chrom, start, end, strand)
#         self.vertex = int(vertex)
#         self.traverse_route = []  # [(HOMOSEQ_REGION, edge_type), ...]

#     def __hash__(self):
#         return hash(self.data)

#     def __repr__(self):
#         return f"{self.chrom}, {self.start}, {self.end}, {self.strand}, {self.size}"

#     def __eq__(self, other):
#         return self.data == other.data

#     def __iter__(self):
#         yield from (self.chrom, self.start + self.ups_rela_start, self.start + self.ups_rela_end, self.strand)

#     def __getitem__(self, index):
#         return self.data[index]

#     def __len__(self):
#         return len(self.data)

#     def fix_coord(self):
#         return (self.chrom, self.start + self.ups_rela_start, self.start + self.ups_rela_end, self.strand)

#     def qnode_relative_region(self, qnode_tuple, logger=logger):
#         """Calculates relative coordinates to a query node (qnode).

#         Args:
#             qnode_tuple: (chrom, start, end, strand) of the query node.
#             logger: Logger object.

#         Returns:
#             (qnode_rela_start, qnode_rela_end) or (np.nan, np.nan) on error.
#         """

#         # Find qnode in traverse_route.
#         qnode_index = -1
#         for i, (node, _) in enumerate(self.traverse_route):
#             if node.data == qnode_tuple:
#                 qnode_index = i
#                 break
#         if qnode_index == -1:
#             logger.error(f"Qnode {qnode_tuple} not in traverse route of {self.data}: {self.traverse_route}")
#             return np.nan, np.nan

#         # Iterate backwards from the qnode.
#         current_node = self
#         qnode_rela_start, qnode_rela_end = self.rela_start, self.rela_end

#         for upstream_node, edge_type in reversed(self.traverse_route[qnode_index:]):
#             if edge_type == "segmental_duplication":
#                 if upstream_node.strand == current_node.strand:
#                     qnode_rela_start = max(0, qnode_rela_start)
#                     qnode_rela_end = min(qnode_rela_end, upstream_node.size)
#                 else:
#                     qsize = qnode_rela_end - qnode_rela_start
#                     qnode_rela_end = min(current_node.size - qnode_rela_start, upstream_node.size)
#                     qnode_rela_start = max(0, qnode_rela_end - qsize)
#             else: #Adjacent
#                 qnode_abs_start = qnode_rela_start + current_node.start
#                 qnode_abs_end = qnode_rela_end + current_node.start
#                 qnode_rela_start = max(0, qnode_abs_start - upstream_node.start)
#                 qnode_rela_end = min(qnode_rela_start + (qnode_abs_end - qnode_abs_start), upstream_node.size)

#             if qnode_rela_end <= qnode_rela_start:
#                 logger.warning(f"Invalid relative coords: start={qnode_rela_start}, end={qnode_rela_end}. "
#                                f"Current: {current_node}, Upstream: {upstream_node}, Edge: {edge_type}")
#                 return np.nan, np.nan

#             current_node = upstream_node

#         return qnode_rela_start, qnode_rela_end