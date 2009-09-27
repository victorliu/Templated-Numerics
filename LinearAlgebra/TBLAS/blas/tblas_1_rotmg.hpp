template <typename ScalarRealT>
void TBLAS_NAME(rotmg,ROTMG)(
  ScalarRealT &d1,
  ScalarRealT &d2,
  ScalarRealT &b1,
  const ScalarRealT &b2,
  ScalarRealT P[5]){
	
	ScalarRealT D1 = d1, D2 = d2, x = b1, y = b2;
	ScalarRealT h11,h12,h21,h22, u, c, s;
	const ScalarRealT G = 4096, G2 = G*G;
	
	
	// Taken from GSL
	  /* case of d1 < 0, appendix A, second to last paragraph */

  if (D1 < 0.0) {
    P[0] = -1;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
    P[4] = 0;
    d1 = 0;
    d2 = 0;
    b1 = 0;
    return;
  }

  if (D2 * y == 0.0) {
    P[0] = -2;                  /* case of H = I, page 315 */
    return;
  }

  c = fabs(D1 * x * x);
  s = fabs(D2 * y * y);

  if (c > s) {
    /* case of equation A6 */

    P[0] = 0.0;

    h11 = 1;
    h12 = (D2 * y) / (D1 * x);
    h21 = -y / x;
    h22 = 1;

    u = 1 - h21 * h12;

    if (u <= 0.0) {             /* the case u <= 0 is rejected */
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }

    D1 /= u;
    D2 /= u;
    x *= u;
  } else {
    /* case of equation A7 */

    if (D2 * y * y < 0.0) {
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }

    P[0] = 1;

    h11 = (D1 * x) / (D2 * y);
    h12 = 1;
    h21 = -1;
    h22 = x / y;

    u = 1 + h11 * h22;

    D1 /= u;
    D2 /= u;

    {
		TBLAS_SWAP(D2,D1);
    }

    x = y * u;
  }

  /* rescale D1 to range [1/G2,G2] */

  while (D1 <= 1.0 / G2 && D1 != 0.0) {
    P[0] = -1;
    D1 *= G2;
    x /= G;
    h11 /= G;
    h12 /= G;
  }

  while (D1 >= G2) {
    P[0] = -1;
    D1 /= G2;
    x *= G;
    h11 *= G;
    h12 *= G;
  }

  /* rescale D2 to range [1/G2,G2] */

  while (fabs(D2) <= 1.0 / G2 && D2 != 0.0) {
    P[0] = -1;
    D2 *= G2;
    h21 /= G;
    h22 /= G;
  }

  while (fabs(D2) >= G2) {
    P[0] = -1;
    D2 /= G2;
    h21 *= G;
    h22 *= G;
  }

  *d1 = D1;
  *d2 = D2;
  *b1 = x;

  if (P[0] == -1.0) {
    P[1] = h11;
    P[2] = h21;
    P[3] = h12;
    P[4] = h22;
  } else if (P[0] == 0.0) {
    P[2] = h21;
    P[3] = h12;
  } else if (P[0] == 1.0) {
    P[1] = h11;
    P[4] = h22;
  }


}
