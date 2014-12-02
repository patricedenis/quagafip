#ifndef _MINMAX_H_
#define _MINMAX_H_

#ifndef quagafip_MIN
#define quagafip_MIN(a,b) ((a)<(b)) ? (a) : (b)
#endif

#ifndef quagafip_MAX
#define quagafip_MAX(a,b) ((a)>(b)) ? (a) : (b)
#endif

#ifndef quagafip_Min
#define quagafip_Min(a,b) ((a)<(b)) ? (a) : (b)
#endif

#ifndef quagafip_Max
#define quagafip_Max(a,b) ((a)>(b)) ? (a) : (b)
#endif

#ifndef quagafip_MIN3
#define quagafip_MIN3(a,b,c) ((quagafip_MIN(a,b)) < (c)) ? (quagafip_MIN(a,b)) : (c)
#endif

#ifndef quagafip_MAX3
#define quagafip_MAX3(a,b,c) ((quagafip_MAX(a,b)) > (c)) ? (quagafip_MAX(a,b)) : (c)
#endif

#endif
