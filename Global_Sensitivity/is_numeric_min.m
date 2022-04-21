% Check a value for being numeric and being >= a specified minimum value.
function  is_numeric_min(value, min, valueName, ident)


if ~isnumeric(value)
    msg = sprintf('%s is not numeric', valueName);
    exception = MException(ident, msg);
    throw(exception);
end

if value < min
    msg = sprintf('%s value of %d is less than %d', valueName, value, min);
    exception = MException(ident, msg);
    throw(exception);
end

end
