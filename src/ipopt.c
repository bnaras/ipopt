#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdbool.h>
#include <string.h>
#define IPOPT_INTTYPES_MACROS 1

#include <IpStdCInterface.h>

typedef struct {
  SEXP eval_f;
  SEXP eval_grad_f;
  SEXP eval_g;
  SEXP eval_jac_g;
  SEXP eval_h;
  SEXP rho;
  SEXP jacobian_structure;
  SEXP hessian_structure;
} RIpoptUserData;

static SEXP make_numeric_from_ipopt(ipindex n, ipnumber *x) {
  SEXP out = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(out), x, sizeof(double) * (size_t) n);
  UNPROTECT(1);
  return out;
}

static SEXP call_r_numeric(SEXP fun, SEXP rho, SEXP x) {
  SEXP call = PROTECT(lang2(fun, x));
  SEXP result = PROTECT(eval(call, rho));
  SEXP numeric_result = PROTECT(coerceVector(result, REALSXP));
  UNPROTECT(3);
  return numeric_result;
}

static bool r_eval_f(ipindex n, ipnumber *x, bool new_x,
                     ipnumber *obj_value, UserDataPtr user_data) {
  (void) new_x;
  RIpoptUserData *ud = (RIpoptUserData *) user_data;
  SEXP rx = PROTECT(make_numeric_from_ipopt(n, x));
  SEXP val = PROTECT(call_r_numeric(ud->eval_f, ud->rho, rx));
  bool ok = LENGTH(val) >= 1;
  if (ok) *obj_value = REAL(val)[0];
  UNPROTECT(2);
  return ok;
}

static bool r_eval_grad_f(ipindex n, ipnumber *x, bool new_x,
                          ipnumber *grad_f, UserDataPtr user_data) {
  (void) new_x;
  RIpoptUserData *ud = (RIpoptUserData *) user_data;
  SEXP rx = PROTECT(make_numeric_from_ipopt(n, x));
  SEXP val = PROTECT(call_r_numeric(ud->eval_grad_f, ud->rho, rx));
  bool ok = LENGTH(val) == n;
  if (ok && LENGTH(val) == n) {
    for (ipindex i = 0; i < n; ++i) grad_f[i] = REAL(val)[i];
  }
  UNPROTECT(2);
  return ok;
}

static bool r_eval_g(ipindex n, ipnumber *x, bool new_x,
                     ipindex m, ipnumber *g, UserDataPtr user_data) {
  (void) new_x;
  RIpoptUserData *ud = (RIpoptUserData *) user_data;
  SEXP rx = PROTECT(make_numeric_from_ipopt(n, x));
  SEXP val = PROTECT(call_r_numeric(ud->eval_g, ud->rho, rx));
  bool ok = LENGTH(val) == m;
  if (ok && LENGTH(val) == m) {
    for (ipindex i = 0; i < m; ++i) g[i] = REAL(val)[i];
  }
  UNPROTECT(2);
  return ok;
}

static bool r_eval_jac_g(ipindex n, ipnumber *x, bool new_x,
                         ipindex m, ipindex nele_jac,
                         ipindex *iRow, ipindex *jCol, ipnumber *values,
                         UserDataPtr user_data) {
  (void) n;
  (void) new_x;
  (void) m;
  RIpoptUserData *ud = (RIpoptUserData *) user_data;
  SEXP jac_struct = ud->jacobian_structure;
  if (values == NULL) {
    int *p = INTEGER(jac_struct);
    for (ipindex k = 0; k < nele_jac; ++k) {
      iRow[k] = p[k] - 1;
      jCol[k] = p[k + nele_jac] - 1;
    }
    return true;
  }

  SEXP rx = PROTECT(make_numeric_from_ipopt(n, x));
  SEXP val = PROTECT(call_r_numeric(ud->eval_jac_g, ud->rho, rx));
  bool ok = LENGTH(val) == nele_jac;
  if (ok && LENGTH(val) == nele_jac) {
    for (ipindex k = 0; k < nele_jac; ++k) values[k] = REAL(val)[k];
  }
  UNPROTECT(2);
  return ok;
}

static bool r_eval_h(ipindex n, ipnumber *x, bool new_x,
                     ipnumber obj_factor, ipindex m, ipnumber *lambda,
                     bool new_lambda, ipindex nele_hess,
                     ipindex *iRow, ipindex *jCol, ipnumber *values,
                     UserDataPtr user_data) {
  (void) new_x;
  (void) new_lambda;
  RIpoptUserData *ud = (RIpoptUserData *) user_data;
  SEXP hess_struct = ud->hessian_structure;
  if (values == NULL) {
    int *p = INTEGER(hess_struct);
    for (ipindex k = 0; k < nele_hess; ++k) {
      iRow[k] = p[k] - 1;
      jCol[k] = p[k + nele_hess] - 1;
    }
    return true;
  }

  SEXP rx = PROTECT(make_numeric_from_ipopt(n, x));
  SEXP robj = PROTECT(ScalarReal(obj_factor));
  SEXP rlambda = PROTECT(allocVector(REALSXP, m));
  for (ipindex i = 0; i < m; ++i) REAL(rlambda)[i] = lambda[i];
  SEXP call = PROTECT(lang4(ud->eval_h, rx, robj, rlambda));
  SEXP val = PROTECT(eval(call, ud->rho));
  SEXP rval = PROTECT(coerceVector(val, REALSXP));
  bool ok = LENGTH(rval) == nele_hess;
  if (ok) {
    for (ipindex k = 0; k < nele_hess; ++k) values[k] = REAL(rval)[k];
  }
  UNPROTECT(6);
  return ok;
}

static void set_options(IpoptProblem prob, SEXP options) {
  if (options == R_NilValue || LENGTH(options) == 0) return;
  SEXP names = getAttrib(options, R_NamesSymbol);
  for (R_xlen_t i = 0; i < XLENGTH(options); ++i) {
    const char *key = CHAR(STRING_ELT(names, i));
    SEXP val = VECTOR_ELT(options, i);
    if (TYPEOF(val) == STRSXP) {
      AddIpoptStrOption(prob, (char *) key, (char *) CHAR(STRING_ELT(val, 0)));
    } else if (TYPEOF(val) == LGLSXP || TYPEOF(val) == INTSXP) {
      AddIpoptIntOption(prob, (char *) key, asInteger(val));
    } else {
      AddIpoptNumOption(prob, (char *) key, asReal(val));
    }
  }
}

SEXP C_ipopt_version(void) {
  int major = 0, minor = 0, release = 0;
  GetIpoptVersion(&major, &minor, &release);
  SEXP out = PROTECT(allocVector(INTSXP, 3));
  INTEGER(out)[0] = major;
  INTEGER(out)[1] = minor;
  INTEGER(out)[2] = release;
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("major"));
  SET_STRING_ELT(names, 1, mkChar("minor"));
  SET_STRING_ELT(names, 2, mkChar("release"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP C_ipopt_solve(SEXP x0, SEXP lower, SEXP upper,
                   SEXP constraint_lower, SEXP constraint_upper,
                   SEXP eval_f, SEXP eval_grad_f, SEXP eval_g,
                   SEXP eval_jac_g, SEXP jacobian_structure,
                   SEXP eval_h, SEXP hessian_structure, SEXP options) {
  ipindex n = LENGTH(x0);
  ipindex m = LENGTH(constraint_lower);
  ipindex nele_jac = nrows(jacobian_structure);
  ipindex nele_hess = nrows(hessian_structure);

  SEXP rho = PROTECT(R_NewEnv(R_EmptyEnv, 1, 29));

  RIpoptUserData ud = {
    eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h, rho,
    jacobian_structure, hessian_structure
  };

  IpoptProblem prob = CreateIpoptProblem(
    n, REAL(lower), REAL(upper), m, REAL(constraint_lower),
    REAL(constraint_upper), nele_jac, nele_hess, 0,
    r_eval_f, r_eval_g, r_eval_grad_f, r_eval_jac_g,
    nele_hess > 0 ? r_eval_h : NULL
  );
  if (prob == NULL) {
    UNPROTECT(1);
    error("failed to create Ipopt problem");
  }

  set_options(prob, options);

  SEXP x = PROTECT(duplicate(x0));
  SEXP g = PROTECT(allocVector(REALSXP, m));
  SEXP mult_g = PROTECT(allocVector(REALSXP, m));
  SEXP mult_x_L = PROTECT(allocVector(REALSXP, n));
  SEXP mult_x_U = PROTECT(allocVector(REALSXP, n));
  memset(REAL(mult_g), 0, sizeof(double) * (size_t) m);
  memset(REAL(mult_x_L), 0, sizeof(double) * (size_t) n);
  memset(REAL(mult_x_U), 0, sizeof(double) * (size_t) n);
  ipnumber obj_val = NA_REAL;

  enum ApplicationReturnStatus status = IpoptSolve(
    prob, REAL(x), REAL(g), &obj_val, REAL(mult_g),
    REAL(mult_x_L), REAL(mult_x_U), &ud
  );
  FreeIpoptProblem(prob);

  SEXP out = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, x);
  SET_VECTOR_ELT(out, 1, ScalarReal(obj_val));
  SET_VECTOR_ELT(out, 2, g);
  SET_VECTOR_ELT(out, 3, ScalarInteger((int) status));
  SET_VECTOR_ELT(out, 4, mult_g);
  SET_VECTOR_ELT(out, 5, mult_x_L);
  SET_VECTOR_ELT(out, 6, mult_x_U);
  SEXP names = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names, 0, mkChar("solution"));
  SET_STRING_ELT(names, 1, mkChar("objective"));
  SET_STRING_ELT(names, 2, mkChar("constraints"));
  SET_STRING_ELT(names, 3, mkChar("status_code"));
  SET_STRING_ELT(names, 4, mkChar("constraint_multipliers"));
  SET_STRING_ELT(names, 5, mkChar("lower_bound_multipliers"));
  SET_STRING_ELT(names, 6, mkChar("upper_bound_multipliers"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(8);
  return out;
}

static const R_CallMethodDef call_methods[] = {
  {"C_ipopt_version", (DL_FUNC) &C_ipopt_version, 0},
  {"C_ipopt_solve", (DL_FUNC) &C_ipopt_solve, 13},
  {NULL, NULL, 0}
};

void R_init_ipopt(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
