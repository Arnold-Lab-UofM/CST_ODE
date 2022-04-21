function  check_function_handle(value, msg, ident)
% Check a value for being a file handle and that the function it refers to
% exists.

if (~isa(value, 'function_handle'))
    msgfmt = 'The %s is not a function handle';
    msg = sprintf(msgfmt, msg);
    throw(MException(ident, msg));
end

handleInfo = functions(value);
if (isempty(handleInfo.file))
    msgfmt = 'The file for function handle %s does not exist.';
    msg = sprintf(msgfmt, handleInfo.function);
    throw(MException(ident, msg));
end

end

