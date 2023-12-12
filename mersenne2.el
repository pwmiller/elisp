;; Mersenne Twister coefficients
(defvar *w* 32)
(defvar *n* 624)
(defvar *m* 397)
(defvar *r* 31)
(defvar *a* #x9908B0DF)
(defvar *u* 11)
(defvar *d* #xFFFFFFFF)
(defvar *s* 7)
(defvar *b* #x9D2C5680)
(defvar *t* 15)
(defvar *c* #xEFC60000)
(defvar *l* 18)
(defvar *f* 1812433253)

;; State array and index
(defvar *MT* (make-vector *n* 0))
(defvar *index* (+ *n* 1))
(defvar *lower-mask* #x7FFFFFFF) ;; (1 << r) - 1
(defvar *upper-mask* #x80000000) ;; lowest w bits of (not lower_mask)

(defun mersenne-seed (seed)
  "Initialize the generator from a seed."
  (aset *MT* 0 seed)
  (dotimes (i (1- *n*))
    (let ((prev (aref *MT* i)))
      (aset
       *MT* (1+ i)
       (logand
        (+ (* *f* (logxor prev (lsh prev (- *w* 2)))) (1+ i)) #xFFFFFFFF))))
  (setq *index* *n*))

(defun mersenne-extract-number ()
  "Extract a tempered value based on MT[index]."
  (when (>= *index* *n*)
    (mersenne-twist)
    (setq *index* 0))
  (let ((y (aref *MT* *index*)))
    (setq y (logxor y (logand (lsh y (- *u*)) *d*)))
    (setq y (logxor y (logand (lsh y *s*) *b*)))
    (setq y (logxor y (logand (lsh y *t*) *c*)))
    (setq y (logxor y (lsh y (- *l*))))
    (incf *index*)
    (logand y #xFFFFFFFF)))

(defun mersenne-twist ()
  "Generate the next n values from the series x_i."
  (dotimes (i *n*)
    (let* ((x
            (+ (logand (aref *MT* i) *upper-mask*)
               (logand (aref *MT* (mod (1+ i) *n*)) *lower-mask*)))
           (xA (lsh x -1)))
      (when (/= (mod x 2) 0)
        (setq xA (logxor xA *a*)))
      (aset *MT* i (logxor (aref *MT* (mod (+ i *m*) *n*)) xA)))))

(defun mersenne-random ()
  "Return a random number in the range [0, 1)."
  (/ (mersenne-extract-number) (expt 2.0 32)))

(defun mersenne-randint (a b)
  "Return a random integer in the range [a, b)."
  (let ((range (- b a))
        (random-value (mersenne-random)))
    (+ a (floor (* random-value range)))))

(defun mersenne-shuffle (list)
  "Shuffle the given list."
  (let ((new-list (copy-sequence list))
        (i (1- (length list))))
    (while (> i 0)
      (let* ((j (mersenne-randint 0 (1+ i)))
             (temp (nth i new-list)))
        (setf (nth i new-list) (nth j new-list))
        (setf (nth j new-list) temp))
      (setq i (1- i)))
    new-list))

(provide 'mersenne)
