function retval = parse_varargin( x )
%PARSE_VARARGIN Summary of this function goes here
%   Detailed explanation goes here
    retval = [];
    if mod(length(x),2)~=0
        junk = 0;
    end
    assert( mod(length(x),2)==0,'parse_varargin error : Length of variable arguement must be EVEN');
    for idx = 1:2:length(x)
        assert(ischar(x{idx}),'Variable arguements have to be name / value pair: First of pair is not char');
        if ischar(x{idx+1})
            s = ['retval.' x{idx} '=''' x{idx+1} ''';'];
        elseif isnumeric(x{idx+1}) || isa(x{idx+1},'logical')
            s = ['retval.' x{idx} '=' num2str(x{idx+1}) ';'];
        end
        try
            eval(s);
        catch
            junk = 0;
        end
    end
end

