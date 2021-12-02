class PassiveElements:
    def __init__(self, passive_coord, centroids, ind_dofs) -> None:
        
        if passive_coord is not None:
            self.elements = self.get_passive_el(passive_coord, centroids)
            self.ind_passive = ind_dofs[self.elements, :]
        else:
            self.elements = None
            self.ind_passive = None

    def get_passive_el(self, passive_coord, centroids):
        """ Gets index of passive elements .
        
        Args:
            passive_coord (:obj:`tuple`): Region that the shape will not be changed.
            centroids (:obj:`numpy.array`): Coordinate (x,y) of the centroid of each element.
            
        Returns:
            Index of passive elements.
        """
        mask = (centroids[:, 0] >= passive_coord[0][0]) & (centroids[:, 0] <= passive_coord[0][1]) & (centroids[:, 1] >= passive_coord[1][0]) & (centroids[:, 1] <= passive_coord[1][1])
        return (mask > 0).nonzero()[0]