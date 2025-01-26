class OpticalObject:
    def __init__(self, Type, Geometry):
        """
        A simple container for an optical object's Type and Geometry.

        Parameters
        ----------
        Type : Any
            The "Type" of the optical object (e.g., "Lens", "Mirror", etc.).
        Geometry : Any
            The geometry data or object describing how the optical object is shaped.
        """
        self.Type = Type
        self.Geometry = Geometry
