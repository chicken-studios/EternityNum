-- How it works
-- The internal representation is as follows: EN.new(sign,layer,mag) == sign*10^10^10^ ... (layer times) mag. So a layer 0 number is just sign*mag, a layer 1 number is sign*10^mag, a layer 2 number is sign*10^10^mag, and so on. If layer > 0 and mag < 0, then the number's exponent is negative, e.g. sign*10^-10^10^10^ ... mag.
-- sign is -1, 0 or 1.
-- layer is a non-negative integer.
-- mag is a Number, normalized as follows: if it is above expl, log10(mag) it and increment layer. If it is below ldown and layer > 0, 10^mag it and decrement layer. At layer 0, sign is extracted from negative mags. Zeroes (self.Sign == 0 or (self.Exp == 0 and self.Layer == 0)) become 0, 0, 0 in all fields. Any infinities have both mag and layer as positive Infinity.
-- Configuration
local expl = 9007199254740991
local ldown = math.log(expl, 10)
local msd = 100
local allow_over = false
local default_digits = 2

---Uses the Lanczos approximation to calculate the Gamma function.
---@param z number
---@return number
function f_gamma(z)
    local g = 7
    local n = 9
    local p = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    }
    local y
    if z < 0.5 then
        y = math.pi / (math.sin(math.pi*z) * f_gamma(1 - z))
    else
        z = z - 1
        local x = p[1]
        for i = 2, #p do
            x = x + (p[i] / (z + i-1))
        end
        local t = z + g + 0.5
        y = math.sqrt(2*math.pi) * t ^ (z + 0.5) * math.exp(-t) * x
    end
    return y
end

---Uses Newton’s method to calculate the Lambert W function.
---@param z number
---@return number
function f_lambertw(z)
    local e = math.exp(1)
    if z == 0 then
        return 0
    end
    if z < -1/e then
        error("No solution for z < -1/e!")
    end

    local tolerance = 1e-10
    local w, wn = nil, nil
    if z <= 0 then
        w = (e*z*math.log(1+math.sqrt(1+e*z)))/(1+e*z + math.sqrt(1+e*z))
    elseif z < e then
        w = z/e
    else
        w = math.log(z)-math.log(math.log(z))
    end

    for i = 1, 100 do
        local ew = math.exp(w)
        local wew = w*ew
        wn = w - (wew - z)/(wew + ew)
        if math.abs(wn - w) < tolerance * math.max(1, math.abs(wn)) then
            return wn
        else
            w = wn
        end
    end
    error("Iteration failed to converge at "..tostring(z).."! (final value is "..tostring(wn)..")")
end


---@class EN
---@field Sign number
---@field Layer number
---@field Exp number
local EN = {}

---@param Sign number
---@param Layer number
---@param Exp number
---@return EN
function Cnew(Sign, Layer, Exp)
    return {Sign=Sign,Layer=Layer,Exp=Exp}
end

local ZERO = Cnew(0,0,0)
local ONE = Cnew(1,0,1)
local nan = Cnew(1,-1,1)
local inf = Cnew(1,math.huge,math.huge)
local DefReturn = ZERO

---Checks if an EternityNum is nan.
---@param value EN
---@return boolean
function EN.isnan(value)
    return value.Sign == nan.Sign and value.Layer == nan.Layer and value.Exp == nan.Exp
end

---Checks if an EternityNum is not finite.
---@param value EN
---@return boolean
function EN.isinf(value)
    return value.Layer == math.huge or value.Exp == math.huge
end

---Checks if an EternityNum is zero.
---@param value EN
---@return boolean
function EN.isZero(value)
    return value.Sign == 0 or (value.Exp == 0 and value.Layer == 0)
end

---Gets the sign of a number.
---@param value number
---@return integer
function math.sign(value)
    if value == 0 then return 0 end
    return value/math.abs(value)
end

local function escapePattern(text)
    return text:gsub("(%p)", "%%%1")
end

---Splits a string into a table using a separator
---@param str string 
---@param sep string
---@return string[] parts
local function split(str, sep)
    local parts = {}
    sep = escapePattern(sep)

    -- iterate through all chunks that are NOT the separator
    for part in string.gmatch(str, "([^" .. sep .. "]+)") do
        table.insert(parts, part)
    end

    return parts
end

---Errorcorrects an EternityNum.
---@param EtNum EN
---@return EN
function EN.errorcorrect(EtNum,exl)
    exl = exl or expl
    local ldwn = math.log(exl,10)
    if EN.isnan(EtNum) then return nan end
    if EN.isinf(EtNum) then return inf end
    if EN.isZero(EtNum) then return ZERO end
    local Sign = EtNum.Sign
    local Layers = EtNum.Layer
    local Exp = EtNum.Exp
    if Layers == 0 and Exp < 0 then
        Exp =- Exp
        Sign =- Sign
    end
    if Layers == 0 and Exp < 1/exl then
        Layers = Layers + 1
        Exp = math.log(Exp,10)
        return Cnew(Sign,Layers,Exp)
    end
    local absExp = math.abs(Exp)
    local signExp = math.sign(Exp)

    while absExp >= exl do
        Layers = Layers + 1
        Exp = signExp * math.log(absExp,10)
        absExp = math.abs(Exp)
        signExp = math.sign(Exp) 
    end
    if absExp < exl then
        while absExp < ldwn and Layers > 0 do
            Layers = Layers - 1
            if Layers == 0 then
                Exp = 10^Exp
            else
                Exp = signExp*10^absExp
                absExp = math.abs(Exp)
                signExp = math.sign(Exp) 
            end
        end
    end
    return Cnew(Sign, Layers, Exp)
end

---Creates a new EternityNum.
---@param Sign number
---@param Layer number
---@param Exp number
---@return EN
function EN.new(Sign, Layer, Exp)
    return EN.errorcorrect({Sign=Sign,Layer=Layer,Exp=Exp})
end

---Converts a number to EternityNum.
---@param value number
---@return EN
function EN.fromNumber(value)
    return EN.new(math.sign(value),0,math.abs(value))
end

---Converts from XeY to EternityNum.
---@param value string
---@return EN
function EN.fromScientific(value)
    local slice = split(value,"e")
    local man = tonumber(slice[1])
    local exp = tonumber(slice[2])
    local sign = math.sign(man)
    
    -- normalize
    local ovf = math.floor(math.log(math.abs(man),10))
    if ovf > 0 then
        man = man / 10^ovf
        exp = exp + ovf
    end

    if exp == 0 then return EN.new(math.sign(man),0,man) end
    if man == 0 then return ZERO end
    if man < 0 then man =- man end

    if exp < 0 then
        if exp < -100 then return ZERO end
        local exp2 = math.log(man,10)+exp
        return EN.errorcorrect(EN.new(sign,1,exp2))
    end
    local exp2 = math.log(man,10)+exp
    local layers = 1

    if exp2 > expl then
        exp2 = math.log(exp2,10)
        layers = layers + 1
    end

    return EN.errorcorrect(EN.new(sign,layers,exp2))
end

---Converts from X;Y to EternityNum.
---@param value string
---@return EN
function EN.fromDSF(value)
    local slice = split(value,";")
    local sign = math.sign(tonumber(slice[1]))
    if sign == 0 then sign = 1 end
    local layers = math.abs(tonumber(slice[1]))
    local exp = tonumber(slice[2])
    return EN.errorcorrect(EN.new(sign,layers,exp))
end

---Converts from string to EternityNum.
---@param value string
---@return EN
function EN.fromString(value)
    if value:find("e") and not value:find(";") then
        return EN.fromScientific(value)
    elseif value:find(";") then
        return EN.fromDSF(value)
    end
    if value == "nan" then return nan end
    if value == "inf" then return inf end
    if value == "" then return DefReturn end

    return EN.fromNumber(tonumber(value))
end


---Converts EternityNum to string.
---@param value EN
---@return string
function EN.tostring(value)
    if EN.isnan(value) then return "nan" end
    if EN.isinf(value) then return "inf" end
    return value.Sign*value.Layer..";"..value.Exp
end

---Converts any valid type to EternityNum.
---@param input any
---@return EN
function EN.convert(input)
    if type(input) == "number" then
        return EN.fromNumber(input)
    elseif type(input) == "string" then
        return EN.fromString(input)
    elseif type(input) == "table" then
        if #input == 2 then
            return EN.fromScientific(input[1].."e"..input[2])
        elseif #input == 3 then
            return EN.errorcorrect(EN.new(input[1],input[2],input[3]))
        elseif input.Sign then
            return EN.errorcorrect(EN.new(input.Sign,input.Layer,input.Exp))
        end
    end
    warn("Invalid input: returning DefaultReturn...")
    return DefReturn
end

---Tries to convert an EternityNum to a number.
---@param value EN
---@return number
function EN.toNumber(value)
    if value.Layer > 1 then
        if math.sign(value.Exp) == -1 then
            return 0
        end
        return value.Sign * math.huge
    end
    if value.Layer == 0 then
        return value.Sign * value.Exp
    elseif value.Layer == 1 then
        return value.Sign * 10^value.Exp
    end

    return math.log(-1,10)
end

---Takes the absolute value of an EternityNum.
---@param value any
---@return EN
function EN.abs(value)
    value = EN.convert(value)

    if value.Sign == 0 then
        return ZERO
    end
    return EN.new(1,value.Layer,value.Exp)
end

---@param value any
---@param value2 any
---@return EN
function EN.maxAbs(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    if EN.absoluteCompare(value,value2) < 0 then return value2 end
    return value
end

---Negates an EternityNum.
---@param value any
---@return EN
function EN.neg(value)
    value = EN.convert(value)
    return EN.new(-value.Sign, value.Layer, value.Exp)
end

---Compares two absolute EternityNum’s.
---1 = greater than
---0 = equals
----1 = less than
---@param value EN
---@param value2 EN
---@return integer
function EN.absoluteCompare(value,value2)
    local layera
    if value.Exp > 0 then
        layera = value.Layer
    else
        layera =- value.Layer
    end

    local layerb
    if value2.Exp > 0 then
        layerb = value2.Layer
    else
        layerb =- value2.Layer
    end

    if layera > layerb then return 1 end
    if layera < layerb then return -1 end

    if value.Exp > value2.Exp then return 1 end
    if value.Exp < value2.Exp then return -1 end

    return 0
end

---Compares two signed EternityNum’s.
---1 = greater than
---0 = equals
----1 = less than
---@param value EN
---@param value2 EN
---@return integer
function EN.compare(value,value2)
    if value.Sign > value2.Sign then return 1 end
    if value.Sign < value2.Sign then return -1 end
    return value.Sign * EN.absoluteCompare(value,value2)
end

---Checks if EN1 is less than EN2
---@param value any
---@param value2 any
---@return boolean
function EN.lt(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    local c = EN.compare(value,value2)
    return c == -1
end

---Checks if EN1 is greater than EN2
---@param value any
---@param value2 any
---@return boolean
function EN.gt(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    local c = EN.compare(value,value2)
    return c == 1
end

---Checks if EN1 is equal to EN2
---@param value any
---@param value2 any
---@return boolean
function EN.eq(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    local c = EN.compare(value,value2)
    return c == 0
end

---Checks if EN1 is less than or equal to EN2
---@param value any
---@param value2 any
---@return boolean
function EN.le(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    local c = EN.compare(value,value2)
    return not (c == 1)
end

---Checks if EN1 is greater than or equal to EN2
---@param value any
---@param value2 any
---@return boolean
function EN.ge(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    local c = EN.compare(value,value2)
    return not (c == -1)
end

---Reciprocates an EternityNum.
---@param value any
---@return EN
function EN.recp(value)
    value = EN.convert(value)
    if value.Exp == 0 then return nan end
    if value.Layer == 0 then
        return EN.new(value.Sign,0,1/value.Exp)
    end
    return EN.new(value.Sign,value.Layer,-value.Exp)
end

---Calculates the logarithm of a value with an arbitrary base.
---@param value any
---@param base any
---@return EN
function baseLog(value, base)
    value = EN.convert(value)
    base = EN.convert(base)

    if value.Sign <= 0 or base.Sign <= 0 then return nan end

    if EN.isnan(base) or EN.isnan(value) then
        return nan
    end

    if value.Layer == 0 and base.Layer == 0 then
        return EN.new(value.Sign, 0, math.log(value.Exp, base.Exp))
    end

    return EN.div(EN.log10(value),EN.log10(base))
end

---Calculates the logarithm (or natural logarithm if base not specified) of a value with an arbitrary base.
---@param value any
---@param base any
---@return EN
function EN.log(value,base)
    if base then
        return baseLog(value,base)
    end

    value = EN.convert(value)

    if value.Sign <= 0 then return nan end

    if value.Layer == 0 then
        return EN.new(value.Sign, 0, math.log(value.Exp,10) * math.log(10))
    elseif value.Layer == 1 then
        return EN.new(math.sign(value.Exp), 0, math.abs(value.Exp) * math.log(10))
    elseif value.Layer == 2 then
        return EN.new(math.sign(value.Exp), 1, math.abs(value.Exp) + math.log(math.log(10),10))
    end

    return EN.new(math.sign(value.Exp), value.Layer-1, math.abs(value.Exp))
end

---Calculates the common logarithm of a value.
---@param value any
---@return EN
function EN.log10(value)
    value = EN.convert(value)

    if value.Sign <= 0 then return nan end

    if value.Layer > 0 then
        return EN.new(math.sign(value.Exp), value.Layer-1, math.abs(value.Exp))
    end

    return EN.new(value.Sign, 0, math.log(value.Exp,10))
end

---Calculates the exponential of a value.
---@param value any
---@return EN
function EN.exp(value)
    value = EN.convert(value)
    local float = (value.Layer == 0) and math.exp(value.Sign * value.Exp) or math.huge
    if value.Layer == 0 and float ~= math.huge then
        return EN.fromNumber(float)
    elseif value.Layer == 0 then
        return EN.new(1,1,value.Sign * math.log(math.exp(1),10) * value.Exp)
    elseif value.Layer == 1 then
        return EN.new(1,2,value.Sign * math.log(math.log(math.exp(1),10),10) + value.Exp)
    else
        return EN.new(1, value.Layer+1, value.Sign * value.Exp)
    end
end

---Adds two values.
---@param value any
---@param value2 any
---@return EN
function EN.add(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)

    if EN.isinf(value) or EN.isinf(value2) then return inf end
    if EN.isZero(value) then return value2 end
    if EN.isZero(value2) then return value end

    if value.Sign == -value2.Sign and value.Layer == value2.Layer and value.Exp == value2.Exp then
        return ZERO
    end

    local a,b

    if value.Layer >= 2 or value2.Layer >= 2 then
        return EN.maxAbs(value,value2)
    end

    if EN.absoluteCompare(value, value2) > 0 then
        a = value
        b = value2
    else
        a = value2
        b = value
    end

    if a.Layer == 0 and b.Layer == 0 then
        return EN.fromNumber(a.Sign * a.Exp + b.Sign * b.Exp)
    end

    local layera = a.Layer * math.sign(a.Exp)
    local layerb = b.Layer * math.sign(b.Exp)

    if layera - layerb >= 2 then return a end

    if layera == 0 and layerb == -1 then
        if math.abs(b.Exp - math.log(a.Exp,10)) > msd then
            return a
        else
            local magdif = 10 ^ (math.log(a.Exp, 10) - b.Exp)
            local man = b.Sign + a.Sign * magdif
            return EN.new(math.sign(man), 1, b.Exp + math.log(math.abs(man),10))
        end
    end
    if layera == 1 and layerb == 0 then
        if math.abs(a.Exp - math.log(b.Exp,10)) > msd then return a end
        local magdif = 10 ^ (a.Exp - math.log(b.Exp, 10))
        local man = b.Sign + a.Sign * magdif
        return EN.new(math.sign(man), 1, math.log(b.Exp,10) + math.log(math.abs(man),10))
    end

    if math.abs(a.Exp - b.Exp) > msd then return a end
    local magdif = 10 ^ (a.Exp - b.Exp)
    local man = b.Sign + a.Sign * magdif
    return EN.new(math.sign(man), 1, b.Exp + math.log(math.abs(man),10))
end

---Subtracts two values.
---@param value any
---@param value2 any
---@return EN
function EN.sub(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    return EN.add(value, EN.neg(value2))
end

---Converts an EN to bnum string
---@param value EN
---@return string
function EN.bnumtostr(value)
    if value.Layer > 2 then if allow_over then return "" end return "inf" end
    if value.Layer == 2 and 10^value.Exp == math.huge then return "inf" end
    if EN.isZero(value) then return "0e0" end

    if value.Layer == 0 then
        local man = (value.Exp / 10^math.floor(math.log(value.Exp,10))) * value.Sign
        return man .. "e" .. math.floor(math.log(value.Exp,10))
    elseif value.Layer == 1 then
        local man = (10 ^ (value.Exp - math.floor(value.Exp))) * value.Sign
        return man .. "e" .. math.floor(value.Exp)
    end
    local man = (10 ^ (value.Exp - math.floor(value.Exp))) * value.Sign
    return man .. "e" .. math.floor(value.Exp)
end

---Multiplies two values.
---@param value any
---@param value2 any
---@return EN
function EN.mul(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)

    if EN.isinf(value) or EN.isinf(value2) then return inf end
    if EN.isZero(value) or EN.isZero(value2) then return ZERO end

    if value.Layer == value2.Layer and value.Exp == -value2.Exp then
        return EN.new(value.Sign * value2.Sign, 0, 1)
    end

    local a,b
    if (value.Layer > value2.Layer) or (value.Layer == value2.Layer and math.abs(value.Exp) > math.abs(value2.Exp)) then
        a = value
        b = value2
    else
        a = value2
        b = value
    end

    if a.Layer == 0 and b.Layer == 0 then
        return EN.fromNumber(a.Sign * b.Sign * a.Exp * b.Exp)
    end
    if a.Layer >= 3 or (a.Layer - b.Layer >= 2) then
        return EN.new(a.Sign * b.Sign, a.Layer, a.Exp)
    end
    if a.Layer == 1 and b.Layer == 0 then
        return EN.new(a.Sign * b.Sign, a.Layer, a.Exp + math.log(b.Exp,10))
    end
    if a.Layer == 1 and b.Layer == 1 then
        return EN.new(a.Sign * b.Sign, a.Layer, a.Exp + b.Exp)
    end

    if (a.Layer == 2 and b.Layer == 1) or (a.Layer == 2 and b.Layer == 2) then
        local temp = EN.new(math.sign(b.Exp), b.Layer - 1, math.abs(b.Exp))
        local nmag = EN.add(EN.new(math.sign(a.Exp), a.Layer - 1, math.abs(a.Exp)), temp)
        return EN.new(a.Sign * b.Sign, nmag.Layer + 1, nmag.Sign * nmag.Exp)
    end
    return nan
end

---Divides two values.
---@param value any
---@param value2 any
---@return EN
function EN.div(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    if EN.isZero(value2) then error("Division by zero") end
    return EN.mul(value, EN.recp(value2))
end

---Calculates the absolute common logarithm of a value.
---@param value any
---@return EN
function EN.abslog10(value)
    value = EN.convert(value)

    if EN.isZero(value) then return nan end

    if value.Layer > 0 then
        return EN.new(math.sign(value.Exp), value.Layer-1, math.abs(value.Exp))
    end

    return EN.new(1, 0, math.log(math.abs(value.Exp),10))
end

---Calculates 10 to the exponent of a value.
---@param value any
---@return EN
function EN.pow10(value)
    value = EN.convert(value)

    if EN.isinf(value) then return inf end

    if value.Layer == 0 then
        local nmag = 10 ^ (value.Sign * value.Exp)
        if nmag < math.huge and math.abs(nmag) > 0.1 then
            return EN.new(1, 0, nmag)
        else
            value = EN.new(value.Sign, value.Layer + 1, math.log(value.Exp,10))
        end
    end

    if value.Sign > 0 and value.Exp > 0 then
        return EN.new(value.Sign, value.Layer + 1, value.Exp)
    end

    if value.Sign < 0 and value.Exp > 0 then
        return EN.new(-value.Sign, value.Layer + 1, -value.Exp)
    end

    return ONE
end

---Calculates a value to the exponent of a value.
---@param value any
---@param value2 any
---@return EN
function EN.pow(value, value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    if EN.isZero(value) and EN.isZero(value2) then
        error("Cannot take zero to the power of zero")
    end
    if EN.isZero(value) then
        return ZERO
    end

    if value.Sign == 1 and value.Layer == 0 and value.Exp == 1 then
        return ONE
    end
    if EN.isZero(value2) then
        return ONE
    end
    if value2.Sign == 1 and value2.Layer == 0 and value2.Exp == 1 then
        return value
    end

    local calc = EN.pow10(EN.mul(EN.abslog10(value), value2))
    if value.Sign == -1 and EN.toNumber(value2) % 2 == 1 then
        return EN.neg(calc)
    elseif value.Sign == -1 and EN.toNumber(value2) < 1e20 then
        local component = EN.fromNumber(math.cos(EN.toNumber(value2) * math.pi))
        return EN.mul(calc, component)
    end
    return calc
end

-- YAY, WE GOT THROUGH THE HARD PART!!! YAYYYYYYYYYYYYYYYYY

-- Abbreviations

function math.smod(x,y)
     local mod = math.fmod(x,y)
     while mod < 0 do
          mod = mod + y
     end
     return mod
end

local function normalformat(num,pre,comma)
  local function toCommas(num)
    local formatted = tostring(num)
    local k = 1
    while true do
      formatted, k = formatted:gsub("^(-?%d+)(%d%d%d)", "%1,%2")
      if k == 0 then break end
    end
    return formatted
  end
  if comma then
    return toCommas(string.format("%."..pre.."f",num))
  else
    return string.format("%."..pre.."f",num)
  end
end
local function exponentialformat(log,pre,eng)
  if not eng then
    local man = 10^(log - math.floor(log))
    local exp = math.floor(log)
    if tonumber(normalformat(man,pre,false)) == 10 then
      man = 1
      exp = exp + 1
    end
    return normalformat(man,pre,false).."e"..normalformat(exp,0,true)
  else
    local con = log/3
    local man = 1000^(con - math.floor(con))
    local exp = math.floor(con)*3
    if tonumber(normalformat(man,pre,false)) == 1000 then
      man = 1
      exp = exp + 3
    end
    return normalformat(man,pre,false).."e"..normalformat(exp,0,true)
  end
end
local illions = {
  Temp = {"","K","M","B","T","Qa","Qt","Sx","Sp","Oc","No"},
  Ones =  {"","U","D","T","Qa","Qt","Sx","Sp","Oc","No"},
  Tens =  {"","De","Vt","Tg","Qag","Qtg","Sxg","Spg","Ocg","Nog"},
  Hundreds =  {"","Ct","Dct","Tct","Qg","Qq","Sg","St","Ot","Nt"},
  Powers = {"Mn","Mc","Nn","Pc","Fm","At","Zp","Yc","Xn","Vc","Mec"},
  PowersN = {1 ,  2  , 3  , 4  , 5  , 6  , 7  , 8  , 9  , 10 , 11  , 12}
}

local function Digit(il,dig)
  return math.smod(math.floor(il/(10^(dig-1))),10)
end
local function smallabbv(illion)
  return illions.Ones[Digit(illion,1)+1] .. illions.Tens[Digit(illion,2)+1] .. illions.Hundreds[Digit(illion,3)+1]
end
local function largeabbv(illion)
  local il = illion
  if illion < 1000 then
    return smallabbv(illion)
  else
    local y = smallabbv(math.smod(illion,1000))
    for q,z in ipairs(illions.PowersN) do
      if q == #illions.PowersN then break end
      local tmp = math.smod(math.floor(il / 1000^z),math.floor(1000^(illions.PowersN[q+1] - z) + 0.001)) -- small correction
      local hyp = (y ~= "" and tmp > 0) and "-" or ""
      if tmp > 0 then
        if tmp == 1 then
          y = illions.Powers[q] .. hyp .. y
        else
          y = smallabbv(tmp) .. illions.Powers[q] .. hyp .. y
        end
      end
    end
    return y
  end
end
local function ABBV(il)
  if il < 10 then
    return illions.Temp[il + 2]
  else
    return largeabbv(il)
  end
end

---Abbreviate a value without standard notation.
---@param value any
---@return string
function EN.defaultAbbv(value)
    if EN.isnan(value) then return "nan" end
    if EN.isinf(value) then return "inf" end
    value = EN.convert(value)
    if EN.lt(value,0) then
        return "-"..EN.defaultAbbv(EN.neg(value))
    end
    if EN.lt(value,"1e9") then
        return normalformat(EN.toNumber(value), 2, true)
    elseif EN.lt(value,"1e15;1e15") then
        local log = value.Exp
        local lay = value.Layer
        if log >= 1e9 and lay == 0 then
            log = math.log(log,10)
            lay = lay + 1
        end
        while log < 15 and lay >= 2 do
            log = 10^log
            lay = lay - 1
        end
        local rep
        if lay < 6 then
            rep = string.rep("e",lay-1)
        else
            rep = "e^" .. normalformat(lay - 1, 0, true).." "
        end
        return rep..exponentialformat(log, 2)
    else
        return "f" .. EN.defaultAbbv(value.Layer)
    end
end

---Checks if a value is between x and y.
---@param val any
---@param x any
---@param y any
---@return boolean
function EN.between(val,x,y)
    val = EN.convert(val)
    x = EN.convert(x)
    y = EN.convert(y)
    return EN.gt(val,x) and EN.lt(val,y)
end

---Takes Value root to Value2
---@param value any
---@param value2 any
---@return EN
function EN.yroot(value,value2)
    value = EN.convert(value)
    value2 = EN.convert(value2)
    return EN.pow(value, EN.recp(value2))
end

---Takes the square root of Value
---@param value any
---@return EN
function EN.sqrt(value)
    value = EN.convert(value)
    return EN.yroot(value,2)
end

---Approximates ln(Gamma(z)) for positive values.
---@param z any
---@return EN
local function lnGammaLanczos(z)
    z = EN.convert(z)

    -- Lanczos coefficients, g = 7, n = 9
    local p = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    }

    local half = EN.fromNumber(0.5)
    local one  = EN.fromNumber(1)
    local g    = EN.fromNumber(7)
    local ln2piHalf = EN.fromNumber(0.5 * math.log(2 * math.pi))

    -- Lanczos uses x = z - 1
    local x = EN.sub(z, one)

    local sum = EN.fromNumber(p[1])
    for i = 2, #p do
        local denom = EN.add(x, EN.fromNumber(i - 1))
        sum = EN.add(sum, EN.mul(EN.fromNumber(p[i]), EN.recp(denom)))
    end

    local t = EN.add(x, EN.add(g, half))

    local result = ln2piHalf
    result = EN.add(result, EN.mul(EN.add(x, half), EN.log(t)))
    result = EN.sub(result, t)
    result = EN.add(result, EN.log(sum))

    return result
end

---Calculates the Gamma function.
---@param value any
---@return EN
function EN.gamma(value)
    value = EN.convert(value)

    if EN.isnan(value) then return nan end
    if EN.isinf(value) then return inf end
    if value.Sign <= 0 then return nan end
    if value.Layer >= 2 then
        return EN.exp(value)
    end

    if value.Layer == 0 then
        if value.Exp < 141 then
            return EN.fromNumber(f_gamma(value.Exp))
        end
        return EN.exp(lnGammaLanczos(value))
    end

    return EN.exp(lnGammaLanczos(value))
end

---Calculates the factorial of a value using Gamma(value + 1).
---@param value any
---@return EN
function EN.fact(value)
    value = EN.convert(value)
    return EN.gamma(EN.add(value, 1))
end

---Returns a random EternityNum between x and y.
---Uses EN arithmetic for the range math.
---@param x any
---@param y any
---@return EN
function EN.rand(x, y)
    x = EN.convert(x)
    y = EN.convert(y)

    if EN.isnan(x) or EN.isnan(y) then return nan end

    if EN.gt(x, y) then
        x, y = y, x
    end

    if EN.eq(x, y) then
        return x
    end

    local r = EN.fromNumber(math.random())
    return EN.add(x, EN.mul(EN.sub(y, x), r))
end

---Returns a random EternityNum between x and y on a logarithmic scale.
---Both bounds must be positive.
---@param x any
---@param y any
---@return EN
function EN.exporand(x, y)
    x = EN.convert(x)
    y = EN.convert(y)

    if EN.isnan(x) or EN.isnan(y) then return nan end
    if EN.le(x, ZERO) or EN.le(y, ZERO) then return nan end

    if EN.gt(x, y) then
        x, y = y, x
    end

    if EN.eq(x, y) then
        return x
    end

    local r = EN.fromNumber(math.random())
    local lx = EN.log10(x)
    local ly = EN.log10(y)

    return EN.pow10(EN.add(lx, EN.mul(EN.sub(ly, lx), r)))
end

---Shifts a value by a digit count.
---Positive digits multiply by 10^digits.
---@param value any
---@param digits any
---@return EN
function EN.shiftByDigits(value, digits)
    value = EN.convert(value)
    digits = EN.convert(digits)

    return EN.mul(value, EN.pow10(digits))
end

---Calculates the principal branch of the Lambert W function.
---@param value any
---@return EN
function EN.lambertw(value)
    value = EN.convert(value)

    if EN.isnan(value) then return nan end
    if EN.eq(value, ZERO) then return ZERO end

    local nvalue = EN.toNumber(value)
    if nvalue ~= math.huge and nvalue ~= -math.huge then
        return EN.fromNumber(f_lambertw(nvalue))
    end

    if EN.le(value, EN.fromNumber(-1 / math.exp(1))) then
        return nan
    end

    if EN.gt(value, ZERO) then
        local lv = EN.log(value)
        return EN.sub(lv, EN.log(lv))
    end

    return nan
end

---Evaluates a math expression and returns an EternityNum.
---Supports +, -, *, /, ^ and parentheses.
---@param expr string
---@return EN
function EN.evalExpression(expr)
    if type(expr) ~= "string" then
        return nan
    end

    local ZERO_VALUE = ZERO or EN.convert("0")
    local unpack = table.unpack or unpack

    local function isIdentChar(c)
        return c:match("[%a_]") ~= nil
    end

    local function isNumChar(c)
        return c:match("[%d%.;eE]") ~= nil
    end

    local function isOp(c)
        return c:match("[%+%-%*%/%^]") ~= nil
    end

    local tokens = {}
    local buffer = ""
    local mode = nil -- "num" or "id"

    local function flush()
        if buffer == "" then return end
        if mode == "num" then
            table.insert(tokens, { type = "num", value = buffer })
        elseif mode == "id" then
            table.insert(tokens, { type = "id", value = buffer })
        end
        buffer = ""
        mode = nil
    end

    local i = 1
    while i <= #expr do
        local c = expr:sub(i, i)

        if c:match("%s") then
            flush()

        elseif mode == "id" then
            if isIdentChar(c) or c:match("%d") then
                buffer = buffer .. c
            else
                flush()
                -- reprocess this character
                i = i - 1
            end

        elseif mode == "num" then
            local last = buffer:sub(-1)
            local canSignAfterE = (last == "e" or last == "E")

            if isNumChar(c) then
                buffer = buffer .. c
            elseif (c == "+" or c == "-") and canSignAfterE then
                buffer = buffer .. c
            else
                flush()
                -- reprocess this character
                i = i - 1
            end

        else
            if isIdentChar(c) then
                mode = "id"
                buffer = c
            elseif isNumChar(c) then
                mode = "num"
                buffer = c
            elseif isOp(c) or c == "(" or c == ")" or c == "," then
                table.insert(tokens, c)
            end
        end

        i = i + 1
    end
    flush()

    local precedence = {
        ["+"] = 1,
        ["-"] = 1,
        ["*"] = 2,
        ["/"] = 2,
        ["^"] = 3,
        ["u-"] = 4,
    }

    local rightAssoc = {
        ["^"] = true,
        ["u-"] = true,
    }

    local output = {}
    local ops = {}
    local ctxStack = {}
    local lastType = nil

    local function markContent()
        for j = 1, #ctxStack do
            if ctxStack[j].isFunc then
                ctxStack[j].hasContent = true
            end
        end
    end

    for idx = 1, #tokens do
        local tok = tokens[idx]
        local nextTok = tokens[idx + 1]

        if type(tok) == "table" and tok.type == "num" then
            table.insert(output, EN.convert(tok.value))
            markContent()
            lastType = "num"

        elseif type(tok) == "table" and tok.type == "id" then
            if nextTok == "(" then
                table.insert(ops, { type = "func", value = tok.value })
                lastType = "id"
            else
                local v = EN[tok.value]
                if v ~= nil and type(v) ~= "function" then
                    table.insert(output, v)
                    markContent()
                    lastType = "num"
                else
                    return nan
                end
            end

        elseif tok == "(" then
            table.insert(ops, tok)

            local isFunc = false
            if #ops >= 2 then
                local prev = ops[#ops - 1]
                if type(prev) == "table" and prev.type == "func" then
                    isFunc = true
                end
            end

            table.insert(ctxStack, {
                isFunc = isFunc,
                commas = 0,
                hasContent = false
            })

            lastType = "lparen"

        elseif tok == "," then
            while #ops > 0 and ops[#ops] ~= "(" do
                table.insert(output, table.remove(ops))
            end
            if #ops == 0 then return nan end

            local ctx = ctxStack[#ctxStack]
            if not ctx or not ctx.isFunc then
                return nan
            end

            ctx.commas = ctx.commas + 1
            lastType = "comma"

        elseif tok == ")" then
            while #ops > 0 and ops[#ops] ~= "(" do
                table.insert(output, table.remove(ops))
            end
            if #ops == 0 then return nan end

            table.remove(ops) -- pop "("

            local ctx = table.remove(ctxStack)
            if ctx and ctx.isFunc then
                local funcTok = table.remove(ops)
                if not funcTok or type(funcTok) ~= "table" or funcTok.type ~= "func" then
                    return nan
                end

                funcTok.argc = ctx.hasContent and (ctx.commas + 1) or 0
                table.insert(output, funcTok)
                markContent()
            else
                markContent()
            end

            lastType = "rparen"

        elseif isOp(tok) then
            local op = tok

            if op == "-" and (lastType == nil or lastType == "op" or lastType == "lparen" or lastType == "comma") then
                op = "u-"
            elseif op == "+" and (lastType == nil or lastType == "op" or lastType == "lparen" or lastType == "comma") then
                goto continue
            end

            while #ops > 0 do
                local top = ops[#ops]
                if top == "(" then break end
                if type(top) == "table" and top.type == "func" then break end

                local p1 = precedence[op]
                local p2 = precedence[top]
                if not p1 or not p2 then return nan end

                if (not rightAssoc[op] and p1 <= p2) or (rightAssoc[op] and p1 < p2) then
                    table.insert(output, table.remove(ops))
                else
                    break
                end
            end

            table.insert(ops, op)
            lastType = "op"
        end

        ::continue::
    end

    while #ops > 0 do
        local top = table.remove(ops)
        if top == "(" or top == ")" then
            return nan
        end
        table.insert(output, top)
    end

    local stack = {}

    for _, tok in ipairs(output) do
        if type(tok) == "table" and tok.Sign ~= nil then
            table.insert(stack, tok)

        elseif type(tok) == "string" and precedence[tok] then
            if tok == "u-" then
                local a = table.remove(stack)
                if not a then return nan end
                table.insert(stack, EN.sub(ZERO_VALUE, a))
            else
                local b = table.remove(stack)
                local a = table.remove(stack)
                if not a or not b then return nan end

                if tok == "+" then
                    table.insert(stack, EN.add(a, b))
                elseif tok == "-" then
                    table.insert(stack, EN.sub(a, b))
                elseif tok == "*" then
                    table.insert(stack, EN.mul(a, b))
                elseif tok == "/" then
                    table.insert(stack, EN.div(a, b))
                elseif tok == "^" then
                    table.insert(stack, EN.pow(a, b))
                end
            end

        elseif type(tok) == "table" and tok.type == "func" then
            local argc = tok.argc or 0
            local args = {}

            for i = argc, 1, -1 do
                local v = table.remove(stack)
                if not v then return nan end
                args[i] = v
            end

            local fn = EN[tok.value]
            if type(fn) ~= "function" then
                return nan
            end

            local ok, result = pcall(fn, unpack(args, 1, argc))
            if not ok then
                return nan
            end

            table.insert(stack, result)
        end
    end

    return stack[1] or ZERO_VALUE
end

function EN.tet(a, b, payload)
    a = EN.convert(a)
    b = EN.toNumber(EN.convert(b))
    payload = EN.convert(payload or 1)
    if b < 0 then
        return EN.log(EN.tetration(a,b+1,payload), a)
    end
    if b == 0 then
        return payload
    end

    local result = EN.mul(EN.pow(a, b%1), payload)
    for i = 1, b do
        result = EN.pow(a, result)
    end

    return result
end

function EN.slog(a, y)
    a = EN.convert(a)
    y = EN.convert(y)

    -- base cases
    if EN.eq(y, 1) then
        return EN.convert(0)
    end

    if EN.eq(y, a) then
        return EN.convert(1)
    end

    -- if y < a, result is between 0 and 1 (fractional height region)
    local low = 0
    local high = 1

    -- expand range until we bracket the answer
    while EN.gt(EN.tet(a, high), y) == false do
        high = high + 1
    end

    -- binary search for fractional height
    for i = 1, 50 do
        local mid = (low + high) / 2
        local val = EN.tet(a, mid)

        if EN.gt(val, y) then
            high = mid
        else
            low = mid
        end
    end

    return EN.convert((low + high) / 2)
end

function EN.min(a, b)
    if EN.compare(a, b) <= 0 then
        return a
    end
    return b
end

function EN.max(a, b)
    if EN.compare(a, b) >= 0 then
        return a
    end
    return b
end

function EN.avg(a, b)
    return EN.div(EN.add(a, b), EN.fromNumber(2))
end

function EN.GSK(n)
    local n2 = EN.mul(n, n)
    local nf = EN.fact(n)

    local term1 = EN.pow(n2, nf)
    local term2 = EN.pow(n, n)

    return EN.sub(term1, term2)
end

return EN
