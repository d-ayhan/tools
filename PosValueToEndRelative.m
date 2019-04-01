function OUT = PosValueToEndRelative(V, FeatureLength)
% For each feature, Converst position value to relative length from the
% beining or the end, whichever is smallest. Relative values to the
% begining are positive while relative values to the end are negative.
% Inputs are row arrays with same length.

le = V - FeatureLength - 1;
le = [le, V];
[m, in] = min(abs(le),[],2);
m(in == 1,1) = m(in == 1,1).*-1;
OUT = m;
end