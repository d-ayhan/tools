function OUT = DetermineLinkOrientation(mateLink)
% case 1 : t[chr:pos[
% case 2 : t]chr:pos[
% case 3 : ]chr:pos]t
% case 4 : [chr:pos[t

token = regexpi(mateLink, '^(\w*)(\[|\]).*?(\[|\])(\w*)$','tokens');
t1 = token{1,1}{1};
t2 = token{1,1}{4};
p = token{1,1}{2};
bit = isempty(t1) + (strcmp('[',p)).*2;
switch bit
    case 2
        c = 1;
        t = t1(2:end);
    case 0
        c = 2;
        t = t1(2:end);
    case 1
        c = 3;
        t = t2(1:end-1);
    case 3
        c = 4;
        t = t2(1:end-1);
end
OUT = {c, t};
end