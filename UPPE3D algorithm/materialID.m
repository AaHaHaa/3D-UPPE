function ID = materialID(material)
%MATERIALID It assigns ID to each material
%   Because the code needs to store the 2D material data, using a cell
%   array to store all the names in string/charater takes up too much
%   memory. Therefore, this code assigns each material an ID so that,
%   instead of a cell array, a memory-saving integer array is used.
%   The size of each cell element is 128 bytes, whereas it is only 8 bytes
%   for a double-precision number. The number of materials in this code
%   will be only around ten or so, so I use the minimal int8 that costs
%   only 1 byte.

switch material
    case 'silica'
        ID = 1;
    case 'chalcogenide'
        ID = 2;
    case 'ZBLAN'
        ID = 3;
    otherwise % no Raman
        ID = 0;
end

end

