package util

type Progress struct {
	errs chan error
	done chan struct{}
}

func NewProgress(total int) *Progress {
	p := &Progress{make(chan error), make(chan struct{})}
	go func() {
		completed := 0
		errorCount := 0
		for err := range p.errs {
			if err == nil {
				completed += 1
			} else {
				errorCount += 1
				if FlagQuiet {
					Warnf("%s", err)
				} else {
					Warnf("\r%s                                    \n", err)
				}
			}

			ratio := 100.0 * (float64(completed) / float64(total))
			Verbosef("\r%d of %d jobs complete (%0.2f%% done, %d errors)",
				completed, total, ratio, errorCount)
		}
		Verbosef("\n")
		p.done <- struct{}{}
	}()
	return p
}

func (p *Progress) JobDone(err error) {
	if p == nil {
		return
	}
	p.errs <- err
}

func (p *Progress) Close() {
	if p == nil {
		return
	}
	close(p.errs)
	<-p.done
}
