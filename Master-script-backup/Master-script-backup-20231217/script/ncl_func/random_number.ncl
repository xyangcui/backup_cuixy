;------------------------------------------------------------
;; take n random_number by uniform.
;; lb: low band.
;; hb: high band.
;; n: random numbers.
;; Information
;; author: cui xiangyang release: v1
;------------------------------------------------------------

function random_number(lb:integer,hb:integer,n:integer)
local lb,hb,n,rd
begin
  ;random_setallseed(363738674, 918273645)
  rd = floattointeger(random_uniform(lb,hb,n))

  return rd

end
