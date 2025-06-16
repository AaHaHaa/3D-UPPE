function n2 = solid_n2(material)
%SOLID_N2 Load solid's n2 value

switch material
    case 'silica'
        n2 = 2.3e-20; % m^2/W
    case 'sapphire'
        n2 = 3e-20; % m^2/W
    case {'N-SF11','SF11'}
        n2 = 7e-20; % m^2/W
    otherwise % unknown
        n2 = 0;
end

end

