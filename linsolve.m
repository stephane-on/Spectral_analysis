function [C,r] = linsolve(A,B,varargin)
%MATLAB Code Generation Library Function

%   Limitations:
%   1. Only UT and LT cases are optimized.  All other options are
%   equivalent to using MLDIVIDE.
%   2. The option structure must be a constant.
%   3. Arrays of input structures are not supported, only a scalar
%      structure input.

%   Copyright 2009-2010 The MathWorks, Inc.
%#codegen

eml_invariant(nargin >= 2, ...
    eml_message('Coder:MATLAB:minrhs'));
eml_invariant(nargin <= 3, ...
    eml_message('Coder:toolbox:linsolve_2'), ...
    'IfNotConst','Fail');
eml_invariant(isa(A,'float') && isa(B,'float'), ...
    eml_message('Coder:MATLAB:linsolve_inputType'));
eml_invariant(ndims(A) == 2 && ndims(B) == 2, ...
    eml_message('Coder:MATLAB:linsolve_inputDim'));
mA = cast(size(A,1),eml_index_class);
nA = cast(size(A,2),eml_index_class);
mB = cast(size(B,1),eml_index_class);
nB = cast(size(B,2),eml_index_class);
minszA = min(mA,nA);
RREQ = nargout == 2;
ONE = ones(eml_index_class);
CZERO = eml_scalar_eg(A,B);
CONE = CZERO + 1;
if nargin > 2
    eml_invariant(isstruct(varargin{1}), ...
        eml_message('Coder:toolbox:linsolve_5'), ...
        'IfNotConst','Fail');
    eml_invariant(isscalar(varargin{1}), ...
        eml_message('Coder:toolbox:linsolve_6'), ...
        'IfNotConst','Fail');
    eml_invariant(eml_is_const(varargin{1}), ...
        eml_message('Coder:toolbox:linsolve_7'), ...
        'IfNotConst','Fail');
    eml_prefer_const(varargin);
end
parms = struct( ...
    'LT',uint32(0), ...
    'UT',uint32(0), ...
    'UHESS',uint32(0), ...
    'SYM',uint32(0), ...
    'POSDEF',uint32(0), ...
    'RECT',uint32(0), ...
    'TRANSA',uint32(0));
poptions = struct( ...
    'CaseSensitivity',false, ...
    'StructExpand',true, ...
    'PartialMatching',true);
pstruct = eml_parse_parameter_inputs(parms,poptions,varargin{:});
LT = eml_get_parameter_value(pstruct.LT,false,varargin{:});
UT = eml_get_parameter_value(pstruct.UT,false,varargin{:});
UHESS = eml_get_parameter_value(pstruct.UHESS,false,varargin{:});
SYM = eml_get_parameter_value(pstruct.SYM,false,varargin{:});
POSDEF = eml_get_parameter_value(pstruct.POSDEF,false,varargin{:});
RECT = eml_get_parameter_value(pstruct.RECT,false,varargin{:});
TRANSA = eml_get_parameter_value(pstruct.TRANSA,false,varargin{:});
%   LT  UT  UHESS  SYM  POSDEF  RECT  TRANSA
%   ----------------------------------------
%   T   F   F      F    F       T/F   T/F
%   F   T   F      F    F       T/F   T/F
%   F   F   T      F    F       F     T/F
%   F   F   F      T    T/F     F     T/F
%   F   F   F      F    F       T/F   T/F
eml_invariant( ...
    islogical(LT) && isscalar(LT) && ...
    islogical(UT) && isscalar(UT) && ...
    islogical(UHESS) && isscalar(UHESS) && ...
    islogical(SYM) && isscalar(SYM) && ...
    islogical(POSDEF) && isscalar(POSDEF) && ...
    islogical(RECT) && isscalar(RECT) && ...
    islogical(TRANSA) && isscalar(TRANSA), ...
    eml_message('Coder:toolbox:linsolve_8'), ...
    'IfNotConst','Fail');
eml_invariant( ...
    (~LT || ~(UT || UHESS || SYM || POSDEF)) && ...
    (~UT || ~(LT || UHESS || SYM || POSDEF)) && ...
    (~UHESS || ~(LT || UT || SYM || POSDEF || RECT)) && ...
    (~SYM || ~(LT || UT || UHESS || RECT)) && ...
    (~POSDEF || (SYM && ~(LT || UT || UHESS || RECT))), ...
    eml_message('Coder:MATLAB:linsolve_CombinationOfFieldsNotCurrentlySupported'));
eml_invariant((~TRANSA && mA == mB) || (TRANSA && nA == mB), ...
    eml_message('Coder:MATLAB:dimagree'));
if TRANSA
    mC = mA;
else
    mC = nA;
end
if LT
    % RECT is ignored.
    C = coder.nullcopy(eml_expand(CZERO,[mC,nB]));
    for j = 1:nB
        for i = 1:minszA
            C(i,j) = B(i,j);
        end
        for i = eml_index_plus(minszA,1):mC
            C(i,j) = 0;
        end
    end
    if TRANSA
        % inv(U)*C --> C
        C = eml_xtrsm('L','L','C','N',minszA,nB,CONE,A,ONE,mA,C,ONE,mC);
        if RREQ
            r = crude_rcond_triangular(A);
        elseif any_diag_zero(A)
            warn_singular;
        end
    else
        % inv(L)*C --> C
        C = eml_xtrsm('L','L','N','N',minszA,nB,CONE,A,ONE,mA,C,ONE,mC);
        if RREQ
            r = crude_rcond_triangular(A);
        elseif any_diag_zero(A)
            warn_singular;
        end
    end
elseif UT
    % RECT is ignored.
    C = coder.nullcopy(eml_expand(CZERO,[mC,nB]));
    for j = 1:nB
        for i = 1:minszA
            C(i,j) = B(i,j);
        end
        for i = eml_index_plus(minszA,1):mC
            C(i,j) = 0;
        end
    end
    if TRANSA
        % inv(L)*C --> C
        C = eml_xtrsm('L','U','C','N',minszA,nB,CONE,A,ONE,mA,C,ONE,mC);
        if RREQ
            r = crude_rcond_triangular(A);
        elseif any_diag_zero(A)
            warn_singular;
        end
    else
        % inv(U)*C --> C
        C = eml_xtrsm('L','U','N','N',minszA,nB,CONE,A,ONE,mA,C,ONE,mC);
        if RREQ
            r = crude_rcond_triangular(A);
        elseif any_diag_zero(A)
            warn_singular;
        end
    end
elseif SYM
    eml_invariant(mA == nA, ...
        eml_message('Coder:MATLAB:square'));
    if POSDEF
        % TODO:  xPOTRS, xPOCON
        % Symmetrize to triu(A), since we don't use a symmetric solver yet.
        for j = ONE:nA
            for i = eml_index_plus(j,1):mA
                A(i,j) = conj(A(j,i));
            end
        end
        if TRANSA
            if RREQ
                [C,r] = eml_lusolve(A',B,RREQ);
            else
                C = eml_lusolve(A',B,RREQ);
            end
        else
            if RREQ
                [C,r] = eml_lusolve(A,B,RREQ);
            else
                C = eml_lusolve(A,B,RREQ);
            end
        end
    else
        % TODO:  xSYTRS, xSYTRF, xSYCON
        % Symmetrize to tril(A), since we don't use a symmetric solver yet.
        for j = ONE:nA
            for i = eml_index_plus(j,1):mA
                A(j,i) = conj(A(i,j));
            end
        end
        if TRANSA
            if RREQ
                [C,r] = eml_lusolve(A',B,RREQ);
            else
                C = eml_lusolve(A',B,RREQ);
            end
        else
            if RREQ
                [C,r] = eml_lusolve(A,B,RREQ);
            else
                C = eml_lusolve(A,B,RREQ);
            end
        end
    end
elseif UHESS
    % TODO:  Specialized Hessenberg solver.
    eml_invariant(mA == nA, ...
        eml_message('Coder:MATLAB:square'));
    A = triu(A,-1);
    if TRANSA
        C = eml_lusolve(A',B,RREQ);
    else
        C = eml_lusolve(A,B,RREQ);
    end
    if RREQ
        r = crude_rcond_triangular(lu(A));
    end
else
    if RECT || (nargin == 2 && mA ~= nA)
        if TRANSA
            [C,r] = eml_qrsolve(A',B,RREQ);
        else
            [C,r] = eml_qrsolve(A,B,RREQ);
        end
    else
        eml_invariant(mA == nA, ...
            eml_message('Coder:MATLAB:square'));
        if TRANSA
            if RREQ
                [C,r] = eml_lusolve(A',B,RREQ);
            else
                C = eml_lusolve(A',B,RREQ);
            end
        else
            if RREQ
                [C,r] = eml_lusolve(A,B,RREQ);
            else
                C = eml_lusolve(A,B,RREQ);
            end
        end
    end
end

%--------------------------------------------------------------------------

function p = any_diag_zero(A)
for k = ones(eml_index_class):min(size(A))
    if A(k,k) == zeros(class(A))
        p = true;
        return
    end
end
p = false;

%--------------------------------------------------------------------------

function warn_singular
eml_warning('Coder:MATLAB:singularMatrix');

%--------------------------------------------------------------------------

function r = crude_rcond_triangular(A)
if isempty(A)
    r = eml_guarded_inf(class(A));
else
    mx = abs(A(1));
    mn = mx;
    for k = cast(2,eml_index_class):min(size(A))
        absAkk = abs(A(k,k));
        if absAkk > mx || isnan(absAkk);
            mx = absAkk;
        elseif absAkk < mn
            mn = absAkk;
        end
    end
    r = mn / mx;
end
%--------------------------------------------------------------------------
