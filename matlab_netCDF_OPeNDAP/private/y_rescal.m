function [rescale_var, rescale_att] = y_rescal(rescale_opts)
% Y_RESCAL returns the scalars rescale_var and rescale_att
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, Wed Feb  9 11:20:27 EST 1994
%     Revision $Revision: 1.5 $
%
% DESCRIPTION:
% y_rescal returns the scalars rescale_var and rescale_att.  Returning
% them this way ensures that getnc.m and getnc_s.m have the same
% values of these variables.  Only alter these if you are sure that
% you know what you are doing.
% 
% INPUT:
% rescale_opts: rescale_var and rescale_att are set to the values
%               rescale_opts(1) and rescale_opts(2) respectively.
%
% OUTPUT:
% rescale_var: If this == 1 then a variable read in by getnc.m and
%              getnc_s.m will be rescaled by 'scale_factor' and
%              'add_offset' if these are attributes of the variable.
%              If == 0 then rescaling will not be done.
% rescale_att: If this == 1 then the attributes '_FillValue',
%              'valid_range', 'valid_min' and 'valid_max' read in by
%              getnc.m and getnc_s.m will be rescaled by
%              'scale_factor' and 'add_offset' when applied to the
%              relevant variable.
%              If == 0 then rescaling will not be done.
%
% EXAMPLE:
% Simply type y_rescal at the matlab prompt.
%
% This function is called by: getnc.m, getnc_s.m
% This function calls: None
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: y_rescal.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

if nargin == 0
  rescale_var = 1;
  rescale_att = 1;
elseif nargin == 1  
  rescale_var = rescale_opts(1);
  rescale_att = rescale_opts(2);
else
  s = [ ' number of input arguments = ' int2str(nargin) ];
  error(s)
end
